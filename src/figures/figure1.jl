""" Figure 1: time-course simulation and data. """

concs, _, g1s1, g2s1 = load(189, 1);
_, _, g1s2, g2s2 = load(189, 2);
_, _, g1s3, g2s3 = load(189, 3);

g1m = (g1s1 .+ g1s2 .+ g1s3) ./ 3;
g2m = (g2s1 .+ g2s2 .+ g2s3) ./ 3;

p = [44.184, 1.24076, 0.0692788, 0.0460918, 0.3822, 0.854034, 0.605391, 0.771326, 0.0138293, 0.00183699, 0.000293753, 0.0127534, 0.00011816, 0.0142541, 60.6069, 0.899573, 1.99993, 0.0748216, 1.99326, 0.468332, 1.99864, 1.22536, 0.000141615, 0.0318616, 0.000216899, 8.80158e-7, 0.598489, 0.00110572, 6.68492, 2.05974, 1.99936, 0.167588, 0.507586, 0.316074, 0.248084, 0.826596, 1.6164e-5, 3.10987e-6, 3.55996e-5, 7.73526e-6, 0.0774056, 8.26708e-5, 3.34656, 2.83739, 0.0907361, 0.108245, 1.9758, 1.96985, 1.9993, 0.210137, 0.0690636, 1.30442e-5, 0.0767181, 0.00991078, 6.87891e-5, 1.45086e-5, 18.2253, 1.1841, 1.00505, 0.0735852, 1.97326, 0.783828, 0.45769, 1.99355, 0.0519941, 0.000533671, 0.00204743, 9.52975e-5, 5.23806e-5, 0.0677505, 0.339953, 0.403341, 0.802518, 0.470576, 1.298, 0.423103];
efcs = getODEparams(p, concs);

G1 = zeros(189, 7, 5)
G2 = zeros(189, 7, 5)


t = LinRange(0.0, 95.0, 189)
for k=1:5 # drug number
    for i = 1:7 # concentration number
        G1[:, i, k], G2[:, i, k], _ = predict(efcs[:, i, k], efcs[:, 1, k], t)
    end
end

p1 = DrugResponseModel.plotavg(G1[:, :, 1], G2[:, :, 1], g1m[:, :, 1], g2m[:, :, 1], 1, :false, concs[1, 1])
pc1 =DrugResponseModel.plotavg(G1[:, :, 1], G2[:, :, 1], g1m[:, :, 1], g2m[:, :, 1], 3, :false, concs[3, 1])
pc2 =DrugResponseModel.plotavg(G1[:, :, 1], G2[:, :, 1], g1m[:, :, 1], g2m[:, :, 1], 5, :false, concs[5, 1])
pc3 =DrugResponseModel.plotavg(G1[:, :, 1], G2[:, :, 1], g1m[:, :, 1], g2m[:, :, 1], 7, :false, concs[7, 1])

p2 = DrugResponseModel.plotavg(G1[:, :, 3], G2[:, :, 3], g1m[:, :, 3], g2m[:, :, 3], 1, :false, concs[1, 3])
pg1 =DrugResponseModel.plotavg(G1[:, :, 3], G2[:, :, 3], g1m[:, :, 3], g2m[:, :, 3], 3, :false, concs[3, 3])
pg2 =DrugResponseModel.plotavg(G1[:, :, 3], G2[:, :, 3], g1m[:, :, 3], g2m[:, :, 3], 5, :false, concs[5, 3])
pg3 =DrugResponseModel.plotavg(G1[:, :, 3], G2[:, :, 3], g1m[:, :, 3], g2m[:, :, 3], 7, :false, concs[7, 3])


function figure1()

    setGadflyTheme()

    
    draw(
        SVG("figure1.svg", 12inch, 8inch),
        plotGrid((2, 4), [p1, pc1, pc2, pc3, p2, pg1, pg2, pg3]),
    )
end
