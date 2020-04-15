using LinearAlgebra: pinv
using DSP: conv


"""
	savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer; boundary_mode = :interpolation)
Apply Savitzky-Golay polynomial smoothing to input data `y` using a polynomial of order `polynomial_order` fit to a moving window `window_size` points wide.
 
# References
1. Savitzky, A., & Golay, M. J. E. (1964). Analytical Chemistry, 36(8), 1627–1639. https://doi.org/10.1021/ac60214a047
2. Steinier, J., Termonia, Y., & Deltour, J. (1972). Analytical Chemistry, 44(11), 1906–1909. https://doi.org/10.1021/ac60319a045
3. Press, W. H., & Teukolsky, S. A. (1990). Computers in Physics, 4(6), 669. https://doi.org/10.1063/1.4822961
"""
function savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer)
	# input validity checks
   	@assert isodd(window_size) "Window size must be an odd integer, i.e. fitting 2m + 1 points around the current value."
	@assert polynomial_order < window_size "Polynomial order must be less than the window size."

	# window size is 2m + 1 points
   	m = (window_size - 1) ÷ 2

	# build the Vandermonde design matrix A. Each row corresponds to a point in the fitting window -m:m
	# and each columns correspond to powers in the range 0:polynomial_order
   	fitting_points = -m:m
   	A = Matrix{Float64}(undef, window_size, polynomial_order + 1)
   	for i in 1:window_size, j in 1:polynomial_order + 1
        A[i,j] = fitting_points[i]^(j - 1)
    end

	# for interpolation we'll want the full pseudo-inverse so we can calculate all the fit values at the edges
	# Ap = y
	C = pinv(A)

	# the filter coefficients are the rows of `C`
	filter_coeffs = C[1,:] * factorial(0)

	# convolve with the filter coefficients with a couple extra steps:
	# 1. because of convolution will reverse coefficients we flip before
	# 2. c = conv(a,b) will return a vector of length(c) = length(a) + length(b) - 1 so we chop off the first and last m points
	smoothed = conv(reverse(filter_coeffs), y)[m+1:end-m]

	# for interpolation edge handling calculate the full fits
	# if we are just smoothing then we can use the design and coefficient matrix as is
	AC = A*C
	smoothed[1:m] = (AC*y[1:window_size])[1:m]
	smoothed[end-m+1:end] = (AC*y[end-window_size+1:end])[end-m+1:end]

	return smoothed
end