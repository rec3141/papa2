"""Error model estimation for DADA2."""

import numpy as np

# Transition row names (matching R's ordering)
TRANS_NAMES = [
    "A2A", "A2C", "A2G", "A2T",
    "C2A", "C2C", "C2G", "C2T",
    "G2A", "G2C", "G2G", "G2T",
    "T2A", "T2C", "T2G", "T2T",
]

# Indices of non-self transitions (12 total)
# Order: A2C,A2G,A2T, C2A,C2G,C2T, G2A,G2C,G2T, T2A,T2C,T2G
_NON_SELF = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14]
_SELF = [0, 5, 10, 15]  # A2A, C2C, G2G, T2T

# Mapping: for each base (0-3), which rows are its transitions (in order A,C,G,T)
_BASE_ROWS = {0: [0, 1, 2, 3], 1: [4, 5, 6, 7], 2: [8, 9, 10, 11], 3: [12, 13, 14, 15]}


def _r_nf(n, span):
    return max(1, min(n, int(np.floor(n * span + 1e-5))))


def _design_matrix_1d(dx, degree):
    cols = [np.ones_like(dx)]
    if degree >= 1:
        cols.append(dx)
    if degree >= 2:
        cols.append(dx ** 2)
    return np.column_stack(cols)


def _loess_fit_coeffs_1d(x, y, weights, x_eval, span=0.75, degree=2):
    """Translate the core R loess local fit for the 1D gaussian case.

    This mirrors the main numerical steps of `ehg127` for our use case:
    1. squared-distance neighbor search
    2. tricube neighborhood weights using `rho = dist[nf] * max(1, f)`
    3. weighted centered polynomial design
    4. column equilibration
    5. QR + SVD pseudoinverse solve
    """
    n = len(x)
    nf = _r_nf(n, span)
    dist = (x - x_eval) ** 2
    order = np.lexsort((np.arange(n, dtype=np.int64), dist))
    psi = order[:nf]

    rho = dist[psi[-1]] * max(1.0, span)
    if rho <= 0.0:
        coeffs = np.zeros(degree + 1, dtype=np.float64)
        coeffs[0] = y[psi[0]]
        return coeffs

    normalized = np.sqrt(dist[psi] / rho)
    kernel = np.zeros(nf, dtype=np.float64)
    inside = normalized < 1.0
    kernel[inside] = (1.0 - normalized[inside] ** 3) ** 3
    w = np.sqrt(weights[psi] * kernel)
    if np.all(w == 0.0):
        coeffs = np.zeros(degree + 1, dtype=np.float64)
        coeffs[0] = y[psi[0]]
        return coeffs

    dx = x[psi] - x_eval
    b = _design_matrix_1d(dx, degree) * w[:, None]
    eta = y[psi] * w

    col_norm = np.linalg.norm(b, axis=0)
    col_norm[col_norm == 0.0] = 1.0
    b_scaled = b / col_norm

    q_mat, r_mat = np.linalg.qr(b_scaled, mode="reduced")
    qty = q_mat.T @ eta

    try:
        u_mat, sigma, vt_mat = np.linalg.svd(r_mat, full_matrices=False)
    except np.linalg.LinAlgError:
        coeffs = np.zeros(degree + 1, dtype=np.float64)
        coeffs[0] = np.average(y[psi], weights=np.maximum(w, 1e-300))
        return coeffs

    tol = sigma[0] * (100.0 * np.finfo(np.float64).eps) if sigma.size else 0.0
    dgamma = np.zeros_like(sigma)
    keep = sigma > tol
    if np.any(keep):
        dgamma[keep] = (u_mat[:, keep].T @ qty) / sigma[keep]

    beta_scaled = vt_mat.T @ dgamma
    return beta_scaled / col_norm


def _loess_fit_point_1d(x, y, weights, x_eval, span=0.75, degree=2):
    coeffs = _loess_fit_coeffs_1d(x, y, weights, x_eval, span=span, degree=degree)
    return coeffs[0]


def _interp_margin(x):
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    mu = 0.005 * max(xmax - xmin, 1e-10 * max(abs(xmin), abs(xmax)) + 1e-30)
    return xmin - mu, xmax + mu


def _build_loess_tree_1d(x, span, cell=0.2):
    """Build the 1D interpolation tree used by R's lowesb/lowese path."""
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    n = len(xs)
    fc = max(1, int(np.floor(n * span * cell)))
    left_bound, right_bound = _interp_margin(xs)
    vertex_map = {left_bound: 0, right_bound: 1}
    vertices = [left_bound, right_bound]

    def add_vertex(value):
        value = float(value)
        found = vertex_map.get(value)
        if found is not None:
            return found
        idx = len(vertices)
        vertices.append(value)
        vertex_map[value] = idx
        return idx

    def build(l, u, left_v, right_v):
        if (u - l + 1) <= fc:
            return {"split": None, "left_v": left_v, "right_v": right_v}
        m = (l + u) // 2
        split = float(xs[m])
        if split == vertices[left_v] or split == vertices[right_v]:
            return {"split": None, "left_v": left_v, "right_v": right_v}
        split_v = add_vertex(split)
        return {
            "split": split,
            "left": build(l, m, left_v, split_v),
            "right": build(m + 1, u, split_v, right_v),
        }

    tree = build(0, n - 1, 0, 1)
    return order, xs, np.asarray(vertices, dtype=np.float64), tree


def _collect_leaf_vertices_1d(node, out):
    if node["split"] is None:
        out.add(node["left_v"])
        out.add(node["right_v"])
        return
    _collect_leaf_vertices_1d(node["left"], out)
    _collect_leaf_vertices_1d(node["right"], out)


def _locate_leaf_1d(node, x_eval):
    while node["split"] is not None:
        if x_eval <= node["split"]:
            node = node["left"]
        else:
            node = node["right"]
    return node


def _hermite_interp_1d(x_eval, x0, x1, val0, val1, deriv0, deriv1):
    width = x1 - x0
    if width == 0.0:
        return val0
    h = (x_eval - x0) / width
    phi0 = (1.0 - h) ** 2 * (1.0 + 2.0 * h)
    phi1 = h ** 2 * (3.0 - 2.0 * h)
    psi0 = h * (1.0 - h) ** 2
    psi1 = h ** 2 * (h - 1.0)
    return phi0 * val0 + phi1 * val1 + (psi0 * deriv0 + psi1 * deriv1) * width


def _lowess_fit(x, y, weights, eval_x=None, span=0.75, degree=2, cell=0.2):
    """1D LOESS for the DADA2 quality-score fit.

    This intentionally follows the R loess point-fit logic for the
    `degree=2`, gaussian, one-predictor case used by `loessErrfun`.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if eval_x is None:
        eval_x = x
    eval_x = np.asarray(eval_x, dtype=np.float64)

    out = np.empty(len(eval_x), dtype=np.float64)
    xmin = x.min()
    xmax = x.max()
    order, xs, vertices, tree = _build_loess_tree_1d(x, span=span, cell=cell)
    ys = y[order]
    ws = weights[order]
    used_vertices = set()
    _collect_leaf_vertices_1d(tree, used_vertices)
    vertex_values = {}
    for vidx in used_vertices:
        coeffs = _loess_fit_coeffs_1d(xs, ys, ws, vertices[vidx], span=span, degree=degree)
        deriv = coeffs[1] if len(coeffs) > 1 else 0.0
        vertex_values[vidx] = (coeffs[0], deriv)

    for i, x_eval in enumerate(eval_x):
        if x_eval < xmin or x_eval > xmax:
            out[i] = np.nan
        else:
            leaf = _locate_leaf_1d(tree, float(x_eval))
            left_v = leaf["left_v"]
            right_v = leaf["right_v"]
            x0 = vertices[left_v]
            x1 = vertices[right_v]
            v0, d0 = vertex_values[left_v]
            v1, d1 = vertex_values[right_v]
            out[i] = _hermite_interp_1d(float(x_eval), x0, x1, v0, v1, d0, d1)
    return out


def loess_errfun(trans):
    """Estimate error rates from transition counts using LOESS smoothing.

    Mirrors R's loessErrfun from errorModels.R.

    Args:
        trans: numpy array (16, ncol), transition counts. Rows are transitions
               (A2A,A2C,...,T2T), columns are quality scores.

    Returns:
        numpy array (16, ncol), estimated error rates.
    """
    ncol = trans.shape[1]
    qq = np.arange(ncol, dtype=np.float64)
    est = np.zeros((12, ncol))  # 12 non-self transitions

    MIN_ERROR_RATE = 1e-7
    MAX_ERROR_RATE = 0.25

    idx = 0
    for nti in range(4):  # A, C, G, T
        base_rows = _BASE_ROWS[nti]
        tot = trans[base_rows].sum(axis=0).astype(np.float64)

        for ntj in range(4):
            if nti == ntj:
                continue
            row = nti * 4 + ntj
            errs = trans[row].astype(np.float64)

            # Log10 error rate with pseudocount
            with np.errstate(divide="ignore", invalid="ignore"):
                rlogp = np.log10((errs + 1.0) / tot)
            rlogp[~np.isfinite(rlogp)] = np.nan

            # Find valid (non-NaN) positions
            valid = ~np.isnan(rlogp)
            if valid.sum() < 2:
                est[idx] = MIN_ERROR_RATE
                idx += 1
                continue

            # LOESS fit on valid points
            x_valid = qq[valid]
            y_valid = rlogp[valid]
            w_valid = tot[valid]
            w_valid = np.maximum(w_valid, 1.0)

            pred = _lowess_fit(x_valid, y_valid, w_valid, eval_x=qq, span=0.75)

            # Fill NaN edges by repeating boundary values
            first_valid = np.where(valid)[0][0]
            last_valid = np.where(valid)[0][-1]
            pred[:first_valid] = pred[first_valid]
            pred[last_valid + 1:] = pred[last_valid]

            # Back-transform from log10
            est[idx] = 10.0 ** pred
            idx += 1

    # Clamp
    est = np.clip(est, MIN_ERROR_RATE, MAX_ERROR_RATE)

    # Build full 16-row matrix with self-transitions as complement
    # Row order: A2A,A2C,A2G,A2T, C2A,C2C,C2G,C2T, G2A,G2C,G2G,G2T, T2A,T2C,T2G,T2T
    # Non-self est order: A2C(0),A2G(1),A2T(2), C2A(3),C2G(4),C2T(5), G2A(6),G2C(7),G2T(8), T2A(9),T2C(10),T2G(11)
    err = np.zeros((16, ncol))
    err[1] = est[0]   # A2C
    err[2] = est[1]   # A2G
    err[3] = est[2]   # A2T
    err[0] = 1.0 - est[0:3].sum(axis=0)  # A2A

    err[4] = est[3]   # C2A
    err[6] = est[4]   # C2G
    err[7] = est[5]   # C2T
    err[5] = 1.0 - est[3:6].sum(axis=0)  # C2C

    err[8] = est[6]   # G2A
    err[9] = est[7]   # G2C
    err[11] = est[8]  # G2T
    err[10] = 1.0 - est[6:9].sum(axis=0)  # G2G

    err[12] = est[9]   # T2A
    err[13] = est[10]  # T2C
    err[14] = est[11]  # T2G
    err[15] = 1.0 - est[9:12].sum(axis=0)  # T2T

    return err


def noqual_errfun(trans):
    """Error estimation ignoring quality scores (constant across Q)."""
    ncol = trans.shape[1]
    err = np.zeros((16, ncol))
    for nti in range(4):
        base_rows = _BASE_ROWS[nti]
        tot = trans[base_rows].sum()
        if tot == 0:
            for ntj in range(4):
                err[nti * 4 + ntj] = 0.25
            continue
        for ntj in range(4):
            row = nti * 4 + ntj
            rate = trans[row].sum() / tot
            err[row] = rate
    return err


def inflate_err(err, inflation):
    """Inflate error rates by a factor (prevents premature convergence).

    Mirrors R's inflateErr.
    """
    out = err * inflation / (1.0 + (inflation - 1.0) * err)
    return np.clip(out, 0.0, 1.0)


def get_initial_err(ncol=41):
    """Get a maximum-possible initial error matrix (all 1.0, matching R)."""
    err = np.full((16, ncol), 1.0)
    return err
