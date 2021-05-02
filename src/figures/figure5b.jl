""" Figure 5b to see why gem 10 nM wih palbociclibs doesn't fit well. """

_, popul2, g1s2, g2s2 = load(189, 2)
_, popul3, g1s3, g2s3 = load(189, 3)

# find G1 std and mean ***** data ******
g1S = cat(g1s1, g1s2, g1s3, dims = 4)
g2S = cat(g2s1, g2s2, g2s3, dims = 4)
g1m = mean(g1S, dims = 4) # mean G1
g2m = mean(g2S, dims = 4) # mean G2

GC = JLD.load("GC.jld")["GC"]
GCsim = zeros(2, 189, 6, 5)
gc_singlesim = zeros(2, 189, 5, 3)
t = LinRange(0.0, 96.0, 189)
concs = hcat([0, 25, 50, 100, 250], [0, 5, 10, 17, 30], [0, 25, 50, 100, 250])
x1 = ["control" "gem 10 nM" "gem 10 nM + palbo 25 nM" "gem 10 nM + palbo 50 nM" "gem 10 nM + palbo 100 nM" "gem 10 nM + palbo 250 nM"]
x2 = ["control" "palbo 25 nM" "palbo 50 nM" "palbo 100 nM" "palbo 250 nM"]
x3 = ["control" "gem 5 nM" "gem 10 nM" "gem 17 nM" "gem 30 nM"]
p=[22.6116, 0.906908, 0.076665, 0.00321998, 0.140348, 0.623161, 2.47439, 2.46818, 0.00107404, 0.00202817, 6.81464e-5, 0.00949839, 0.00309843, 0.130674, 1.11947, 0.174963, 2.27961, 0.0117939, 0.00275518, 0.162717, 4.50623, 4.9279, 0.00698243, 0.00333901, 2.85654e-5, 0.0241219, 0.00038035, 0.00244354, 20.2521, 0.20926, 0.449486, 0.00905644, 0.376322, 0.540614, 2.42192, 2.4884, 3.00388e-5, 1.50455e-5, 0.000278125, 0.000434414, 0.0119472, 0.184269, 0.293064, 0.784923, 0.248949, 0.378405, 0.582056, 2.03156]
combinEffects = DrugResponseModel.my_helper(p)
single_effecs = getODEparams(p, concs)

for j=1:5
    for i = 1:6 # concentration number
        GCsim[1, :, i, j], GCsim[2, :, i, j], _ = predict(combinEffects[:, i, j], combinEffects[:, 1, j], t)
    end
end

for i=1:3
    for j=1:5
        gc_singlesim[1, :, j, i], gc_singlesim[2, :, j, i], _ = predict(single_effecs[:, j, i], single_effecs[:, 1, i], t)
    end
end

function plot_separate()

    p1 = DrugResponseModel.plot_fig5(t, GCsim[1, :, :, 3], GC[1, 1:189, :, 3], x1, "Gem 10 + Palbo combo", "G1", "A")
    p2 = DrugResponseModel.plot_fig5(t, GCsim[2, :, :, 3], GC[2, 1:189, :, 3], x1, "Gem 10 + Palbo combo", "G2", "B")
    p3 = DrugResponseModel.plot_fig5(t, gc_singlesim[1, :, :, 3], g1m[:, [1, 4, 5, 6, 7], 5], x2, "Palbo alone", "G1", "C")
    p4 = DrugResponseModel.plot_fig5(t, gc_singlesim[2, :, :, 3], g2m[:, [1, 4, 5, 6, 7], 5], x2, "Palbo alone", "G2", "D")
    p5 = DrugResponseModel.plot_fig5(t, gc_singlesim[1, :, :, 2], g1m[:, [1, 4, 5, 6, 7], 3], x3, "Gem alone", "G1", "E")
    p6 = DrugResponseModel.plot_fig5(t, gc_singlesim[2, :, :, 2], g2m[:, [1, 4, 5, 6, 7], 3], x3, "Gem alone", "G2", "F")
    fig5b = plot(p1, p2, p3, p4, p5, p6, size = (800, 1400), layout = (3, 2))
    savefig(fig5b, "figure5b.svg")
end