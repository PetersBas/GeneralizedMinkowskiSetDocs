---
title: Julia sofware for computing projections onto the generalized Minkowski set
author: |
    Bas Peters^\*^\, Felix J. Herrmann^\#^\
    ^\*^University of British Columbia.\
	^\#^Georgia Institute of Technology
runninghead: SetIntersectionProjection.
classoption:
    - paper
---

##[Preprint paper available online soon]
##[Code on github merged with SetIntersectionProjection](https://github.com/slimgroup/SetIntersectionProjection.jl)

This is the documentation main page corresponding to the **Julia 1.1** software developed by Bas Peters and Felix J. Herrmann that computes projections of vectorized 2D and 3D images/models (``m \in \mathbb{R}^N``) onto a generalization of the Minkowski set:

```math
\mathcal{M} \equiv \{ m = u + v \: | \: u \in \bigcap_{i=1}^p \mathcal{D}_i, \: v \in \bigcap_{j=1}^q \mathcal{E}_j, \: m \in \bigcap_{k=1}^r \mathcal{F}_k \},
```
The vector ``m`` is an element of the generalized Minkowski set if ``m`` is an element of the intersection of ``r`` sets ``\mathcal{F}_k`` and also the sum of two
components ``u \in \mathbb{R}^N`` and ``v \in \mathbb{R}^N``. The vector ``u``
is an element of the intersection of ``p`` sets ``\mathcal{D}_i`` and ``v`` is an
element of the intersection of ``q`` sets ``\mathcal{E}_j``.

This is a set designed to descripe prior knowledge to regularize inverse problems using constraints. Each set describes prior knowledge on the model parameters itself, or, properties of one of the components. See below for examples of sets we can work with. 

An important feature of the algorithms and software, is their focus on problems where each set definition can also include a, possibly non-orthogonal, linear operator. The Euclidean projection onto the generalized Minkowski set may then be formulated using sums of indicator functions as
```math
\min_{u,v,w} \frac{1}{2} \| w - m \|_2^2 + \sum_{i=1}^p \iota_{\mathcal{D}_i}(A_i u) + \sum_{j=1}^q \iota_{\mathcal{E}_j}(B_j v) + \sum_{k=1}^r \iota_{\mathcal{F}_k}(C_k w) + \iota_{w=u+v}(w,u,v),
```
where ``A_i``, ``B_i``, and ``C_i`` are different linear operators of each ``i``, ``j``, and ``k``.

The generalized Minkowski set, algorithms, and software are designed
 
- for applications in imaging inverse problems.
- as a plug-in projector for other algorithms that solve
	``\min_m f(m) \:\: \text{s.t.} \:\: m \in \mathcal{M}``
	, e.g., a (spectral) projected gradient / projected quasi-Newton / projected-Newton method. 
- as a solver for linear inverse problem with a linear forward operator ``F \in \mathbb{R}^{M \times N}``, data constraints such as ``\mathcal{C}^\text{data} = \{ x \: | \: \| Fx - d_\text{observed} \|_2 \leq \sigma \}``, and model property constraints formulated as
	```math
	\min_{x} \frac{1}{2} \| x - m \|_2^2 \quad \text{s.t.} \quad \begin{cases} x \in \mathcal{M} \\
	x \in \mathcal{G}^{\text{data}}
	\end{cases}.
	```

## Applications

 - [Generalized Minkowski decomposition of a video](https://github.com/slimgroup/SetIntersectionProjection.jl/blob/master/examples/GeneralizedMinkowski/Minkowski_video_decomposition.jl)
 - [Generalized Minkowski projection and decompositions of a geological model](https://github.com/slimgroup/SetIntersectionProjection.jl/blob/master/examples/GeneralizedMinkowski/example_2D_Minkowski_projection.jl)
 - [Seismic full-waveform inversion with generalized Minkowski constraints (coming soon)]
 
## Computational features

- this package is based on the [SetIntersecionProjection](https://petersbas.github.io/SetIntersectionProjectionDocs/) package that was designed to compute projections onto intersections of sets. Not all SetIntersecionProjection features are available for generalized Minkowski sets at the moment.
- parametric typing for Float32 and Float64 support
- designed for small 2D up to larger 3D models/grids
- includes scripts to set up projectors and linear operators For 2D and 3D models
- constraints may be defined for the matrix/tensor model and for columns/slices/fibers simultaneously
- linear operators may be: SparseMatrixCSC or [JOLI](https://github.com/slimgroup/JOLI.jl) DCT/DFT/Curvelet matrix-free operators

## List of constraints & linear operators

#### Table: {#set-overview}
|  descriptions | set |
|--- |--- | 
| bounds  | ``\{ m \: | \: l[i] \leq m[i] \leq u[i] \}``| `m_min     = 0.0`,`m_max     = 255.0`, `set_type  = "bounds"`, `TD_OP     = "identity"`, `app_mode  = ("matrix","")`,`custom_TD_OP = ([],false)` |
| transform-domain bounds | ``\{ m \: | \: l[i] \leq (A m)[i] \leq b[i] \}`` | `constraint["use_TD_bounds_1"]=true`, `constraint["TD_LB_1"]=l`, `constraint["TD_UB_1"]=u` `constraint["TDB_operator_1"]=A` |
| (special case) vertical (approximate) monotonicity | ``\{ m \: | \: l[i] \leq (D_z \otimes I_x) m)[i] \leq u[i] \}`` | `constraint["use_TD_bounds_1"]=true`, `constraint["TD_LB_1"]=-eps`, `constraint["TD_UB_1"]=+eps`, `constraint["TDB_operator_1"]=D_z` |
| transform-domain ``\ell_1`` | ``\{ m \: | \: \| A m \|_1 \leq \sigma \}`` | `constraint["use_TD_l1_1"]=true`, `constraint["TD_l1_operator_1"]=A`, `constraint["TD_l1_sigma_1"] = sigma` |
| transform-domain ``\ell_2`` | ``\{ m \: | \: \| A m \|_2 \leq \sigma \}`` | `constraint["use_TD_l2_1"]=true`, `constraint["TD_l2_operator_1"]=A`, `constraint["TD_l2_sigma_1"] = sigma` |
| transform-domain annulus | ``\{ m \: | \: \sigma_l \leq \| A m \|_2 \leq \sigma_u \}`` | `constraint["use_TD_annulus_1"]=true`, `constraint["TD_annulus_operator_1"]=A`, `constraint["TD_annulus_sigma_min_1"] = sigma_l`, `constraint["TD_annulus_sigma_max_1"] = sigma_u` |
| transform-domain cardinality | ``\{ m \: | \: \operatorname{card}(Am) \leq k \}``, ``k`` is a positive integer | `constraint["use_TD_card_1"]=true`, `constraint["TD_card_operator_1``]=A`, `constraint["card_1"]` |
| transform-domain nuclear norm | ``\{ m \: | \: \sum_{j=1}^k \lambda[j] \leq \sigma \}``, with ``Am = \operatorname{vec}( \sum_{j=1}^{k}\lambda[j] u_j v_j^* )`` is the SVD | `constraint["use_TD_nuclear_1_"]=true`, `constraint["TD_nuclear_operator_1"]=A`, `constraint["TD_nuclear_norm_1"]  = sigma` |
| transform-domain rank constraint | ``\{ m \: | \:  Am = \operatorname{vec}( \sum_{j=1}^{r}\lambda[j] u_j v_j^*) \}``, ``r < \text{min}(n_z,n_x)`` | `constraint["use_TD_rank_1"]=true`, `constraint["TD_rank_operator_1"]=A`, `constraint["max_TD_rank_1"]=r` |
| subspace constraints | ``\{ m \: | m = A c, \:\: c \in \mathbb{C}^M \}`` | `constraint["use_subspace"]=true`, `constraint["A"]=A`, `constraint["subspace_orthogonal"]=true`

: Overview of constraint sets that the software currently supports. A new constraint set may be added by providing a projection onto the set (without linear operator) and a sparse linear operator or equivalent matrix-vector product together with its adjoint. Vector entries are indexed as ``m[i]``.


#### Table: {#LOP-overview}
|  descriptions | Linear operator | code
|--- |--- | --- |
|discrete derivative in one direction | ``D_z \otimes I_x`` , ``I_z \otimes D_x`` | "D_z", "D_x" 	|
|discrete derivative in all directions | ``\begin{pmatrix} D_z \otimes I_x \\ I_z \otimes D_x \end{pmatrix}`` | "D2D" or "D3D" | 
| identity matrix | ``I`` | "identity" |
| discrete cosine transform | | "DCT" |
| discrete Fourier transform | | "DFT" |
| curvelet transform | | "curvelet" |
| wavelet transform | | "wavelet" |

: Overview of the linear operators that we currently set up. Software can work with any linear operator as long it is one of the types `SparseMatrixCSC` or `JOLI` operator. Possible conversion to CDS format happens in the software. Operator math is shown for the 2D case. Curvelets require the separate installation of the [CurveLab](http://curvelet.org/software.html) software.

.

```math_def
\def\bb{\mathbf b}
\def\bc{\mathbf c}
\defd{\mathbf d}
\def\bg{\mathbf g}
\def\bh{\mathbf h}
\def\bl{\mathbf l}
\defm{\mathbf m}
\def\bp{\mathbf p}
\def\bq{\mathbf q}
\def\br{\mathbf r}
\def\bs{\mathbf s}
\def\bu{\mathbf u}
\defv{\mathbf v}
\def\bw{\mathbf w}
\defy{\mathbf y}
\defx{\mathbf x}
\def\bz{\mathbf z}
%\def\argmin{\operatornamewithlimits{arg min}}
\def\argmin{\mathop{\rm arg\min}}
```