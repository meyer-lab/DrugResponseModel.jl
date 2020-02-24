""" In this file we fit all the drugs att once. """

""" Calculates the cost function for all of the concentrations. """
function costAll(p, g1, g2, nG1::Int, nG2::Int, nD1::Int, nD2::Int)
    # We are assuming each of the data matrices for G1 and G2 are 192x8x4 for all 4 drugs we have. We could make this size-flexible for later that we want to use it for more drugs.
    t = LinRange(0.0, 95.5, 192)
    G1pred = zeros(192, 8, length(g1[1,1,:]))
    G2pred = zeros(192, 8, length(g2[1,1,:]))
    allNorms = 0.0
    j=1
    for i=1:length(g1[1,1,:])
        G1[:,:,i], G2[:,:,i] = predict(p[j:j+4], g1[1,1,i] + g2[1,1,i], t, nG1, nG2, nD1, nD2)
        tempNorm = norm(G1[:,:,i] - g1[:,:,i]) + norm(G2[:,:,i] - g2[:,:,i])
        allNorms += tempNorm
        j += 5
    end

    return allNorms
end
