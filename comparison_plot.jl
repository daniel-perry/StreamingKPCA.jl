# Random Fourier Features
# from: Rahimi,Recht, "Random Features for Large-Scale Kernel Machines", NIPS 2007
#
# Danny Perry (dperry@cs.utah.edu)
# Oct 2014
#
# compare to Gaussian kernel:

using PyPlot
using PyCall
@pyimport matplotlib.lines as mlines
using Debug
using MAT

using RandomFourierFeatures
using Kernel
using Sketch


macro nonzero(default, list)
	:([$list[find($list .> 0)],$default])
end

X = []
Labels = []

name = "gaussian"

d = 0
n = 0

@debug function showresults(fn,dataset)

	FDstyles = ["bo-","gs-","yD-","mv-","c<-","r+-"]

	rankkstyles = ["bs--","rs--","gs--","ms--","cs--","ys--"]
	Krankkstyles = ["b+:","r+:","g+:","m+:","c+:","y+:"]


	if dataset == "cpu"
		n = 3000
	end

	file = matopen(fn)
	errs_fro = read(file, "errs_fro")
	errs_l2 = read(file, "errs_l2")
	timings = read(file, "timings")
	timings_test = read(file, "timings_test")

	errs_FD_fro = read(file, "errs_FD_fro")
	errs_FD_l2 = read(file, "errs_FD_l2")
	timings_FD = read(file, "timings_FD")
	timings_test_FD = read(file, "timings_test_FD")
	errs_nystrom_fro = read(file, "errs_nystrom_fro")
	errs_nystrom_l2 = read(file, "errs_nystrom_l2")
	timings_nystrom = read(file, "timings_nystrom")
	#timings_test_nystrom = read(file, "timings_test_nystrom")
	errs_nystromwr_fro = read(file, "errs_nystromwr_fro")
	errs_nystromwr_l2 = read(file, "errs_nystromwr_l2")
	timings_nystromwr = read(file, "timings_nystromwr")
	timings_test_nystromwr = read(file, "timings_test_nystromwr")
	errs_rff_nystrom_fro = read(file, "errs_rff_nystrom_fro")
	errs_rff_nystrom_l2 = read(file, "errs_rff_nystrom_l2")
	timings_rff_nystrom = read(file, "timings_rff_nystrom")
	timings_test_rff_nystrom = read(file, "timings_test_rff_nystrom")
	dims = read(file, "dims")
	ells = read(file, "ells")
	epsilons = read(file, "epsilons")
	#sketch_size = read(file, "sketch_size")
  #sketch_size = 5
	#dims = [1:length(errs_fro)]
	dims = dims[1:length(errs_fro)]
	#epsilons = [1e1 1e0 1e-1 1e-2 1e-3]
	#ells = zeros(length(epsilons))
	#for i=1:length(epsilons)
	#	ells[i] = ceil(sketch_size + sketch_size/epsilons[i])
	#end
	close(file)


	marksize = 8
	fontsize = 25
	tickfontsize = 20
	linewidth = 2

	ellis = collect(1:length(ells))
	titlename = dataset
	addlim = 1
	

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims.^2)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims.^2)))
	ymin = min(ymin, minimum(@nonzero(ymin,errs_fro)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_fro)))
	plot(dims.^2,errs_fro,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	#plot(dims.^2,errs_nystrom_fro,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	ymin = min(ymin, minimum(@nonzero(ymin,errs_nystromwr_fro)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_nystromwr_fro)))
	plot(dims.^2,errs_nystromwr_fro,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,ells[li].*dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,ells[li].*dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,errs_FD_fro[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,errs_FD_fro[li,1:length(dims)])))
  	plot(ells[li] .* dims, reshape(errs_FD_fro[li,1:length(dims)],length(dims)) ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.10)
	#legend(names, loc="upper right", frameon = false, numpoints = 1)
  yscale("log")
  xscale("log")
  #ylabel("Kernel Frobenius Error", fontsize = fontsize)
  xlabel("space", fontsize = fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	#title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_space_error_fro.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_space_error_fro.pdf",dataset),bbox_inches="tight")
	

	#figure()
	#l1 = mlines.Line2D([], [], color="r", marker="^", linewidth=linewidth, markersize=marksize, label=names[1], numpoints=1)
	##legend(handles=[l1])
	###legend(names, loc="upper right", frameon = false, numpoints = 1)
	#savefig(@sprintf("#legend.png",dataset),bbox_inches="tight")
	#savefig(@sprintf("#legend.pdf",dataset),bbox_inches="tight")

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims.^2)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims.^2)))
	ymin = min(ymin, minimum(@nonzero(ymin,errs_l2)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_l2)))
	plot(dims.^2,errs_l2 ,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	#plot(dims.^2,errs_nystrom_l2,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	ymin = min(ymin, minimum(@nonzero(ymin,errs_nystromwr_l2)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_nystromwr_l2)))
	plot(dims.^2,errs_nystromwr_l2 ,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,ells[li].*dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,ells[li].*dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,errs_FD_l2[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,errs_FD_l2[li,1:length(dims)])))
  	plot(ells[li] .* dims, reshape(errs_FD_l2[li,1:length(dims)],length(dims)) ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.10)
	legend(names, loc="upper right", frameon = false, numpoints = 1)
  yscale("log")
  xscale("log")
	title(dataset)
  ylabel("Kernel Spectral Error", fontsize=fontsize)
  xlabel("space", fontsize=fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	#ylim(ymin*.90,ymax*1.10)
	title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_space_error_l2.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_space_error_l2.pdf",dataset),bbox_inches="tight")

	



	############ SAMPLE SIZE PLOTS

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims)))
	ymin = min(ymin, minimum(@nonzero(ymin,errs_fro)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_fro)))
	plot(dims,errs_fro ,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	#plot(dims.^2,errs_nystrom_fro,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	ymin = min(ymin, minimum(@nonzero(ymin,errs_nystromwr_fro)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_nystromwr_fro)))
	plot(dims,errs_nystromwr_fro ,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,errs_FD_fro[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,errs_FD_fro[li,1:length(dims)])))
  	plot(dims, reshape(errs_FD_fro[li,1:length(dims)],length(dims)) ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.01)
	legend(names, loc="upper right", frameon = false, numpoints = 1)
  yscale("log")
  xscale("log")
  ylabel("Kernel Frobenius Error", fontsize = fontsize)
  xlabel("sample size", fontsize = fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	#ylim(ymin*.90,ymax*1.10)
	#ylim(1e0,1e2)
	title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_sample_error_fro.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_sample_error_fro.pdf",dataset),bbox_inches="tight")
	

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims)))
	ymin = min(ymin, minimum(@nonzero(ymin,errs_l2)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_l2)))
	plot(dims,errs_l2 ,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	ymin = min(ymin, minimum(@nonzero(ymin,errs_nystromwr_l2)))
	ymax = max(ymax, maximum(@nonzero(ymax,errs_nystromwr_l2)))
	#plot(dims.^2,errs_nystrom_l2,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	plot(dims,errs_nystromwr_l2 ,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,errs_FD_l2[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,errs_FD_l2[li,1:length(dims)])))
  	plot(dims, reshape(errs_FD_l2[li,1:length(dims)],length(dims))  ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.01)
	legend(names, loc="upper right", frameon = false, numpoints = 1)
  yscale("log")
  xscale("log")
	title(dataset)
  ylabel("Kernel Spectral Error", fontsize=fontsize)
  xlabel("sample size", fontsize=fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_sample_error_l2.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_sample_error_l2.pdf",dataset),bbox_inches="tight")
	


	############ RUNTIME PLOTS - training

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims)))
	ymin = min(ymin, minimum(@nonzero(ymin,timings)))
	ymax = max(ymax, maximum(@nonzero(ymax,timings)))
	plot(dims,timings ,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	#plot(dims.^2,errs_nystrom_fro,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	ymin = min(ymin, minimum(@nonzero(ymin,timings_nystromwr)))
	ymax = max(ymax, maximum(@nonzero(ymax,timings_nystromwr)))
	plot(dims,timings_nystromwr, color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,timings_FD[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,timings_FD[li,1:length(dims)])))
  	plot(dims, reshape(timings_FD[li,1:length(dims)],length(dims)) ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.01)
	legend(names, loc="upper right", frameon = false, numpoints = 1)
  yscale("log")
  xscale("log")
  ylabel("Runtime (sec)", fontsize = fontsize)
  xlabel("sample size", fontsize = fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	#ylim(ymin*.90,ymax*1.10)
	#ylim(1e0,1e2)
	title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_timings_train.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_timings_train.pdf",dataset),bbox_inches="tight")
	

	############ RUNTIME PLOTS - tests

	figure()
	xmin=1e6
	xmax=-1e6
	ymin=1e6
	ymax=-1e6
	xmin = min(xmin, minimum(@nonzero(xmin,dims)))
	xmax = max(xmax, maximum(@nonzero(xmax,dims)))
	ymin = min(ymin, minimum(@nonzero(ymin,timings_test)))
	ymax = max(ymax, maximum(@nonzero(ymax,timings_test)))
	plot(dims,timings_test ,color="r",marker="^",linewidth=linewidth,markersize=marksize)
	hold(:on)
	#plot(dims.^2,errs_nystrom_fro,color="k",marker="*",linewidth=linewidth,markersize=marksize)
	ymin = min(ymin, minimum(@nonzero(ymin,timings_test_nystromwr)))
	ymax = max(ymax, maximum(@nonzero(ymax,timings_test_nystromwr)))
	plot(dims,timings_test_nystromwr, color="k",marker="*",linewidth=linewidth,markersize=marksize)
	names = ["RNCA","Nystrom"]
	for li=ellis
		xmin = min(xmin, minimum(@nonzero(xmin,dims)))
		xmax = max(xmax, maximum(@nonzero(xmax,dims)))
		ymin = min(ymin, minimum(@nonzero(ymin,timings_test_FD[li,1:length(dims)])))
		ymax = max(ymax, maximum(@nonzero(ymax,timings_test_FD[li,1:length(dims)])))
  	plot(dims, reshape(timings_test_FD[li,1:length(dims)],length(dims)) ,color=@sprintf("%s",FDstyles[li][1]),marker=@sprintf("%s",FDstyles[li][2]),linewidth=linewidth,markersize=marksize)
		push!(names, @sprintf("SKPCA (%d)",ells[li]))
	end
	ylim(ymin*0.90,ymax*1.10)
	xlim(xmin*0.90,xmax*1.01)
	#legend(names, loc="upper center", frameon = false, numpoints = 1, ncol=length(names))
	legend(names, loc="upper center", frameon = false, numpoints = 1)
  yscale("log")
  #xscale("log")
  ylabel("Runtime (sec)", fontsize = fontsize)
  xlabel("sample size", fontsize = fontsize)
	tick_params(axis="both", which="major", labelsize=tickfontsize)
	tick_params(axis="both", which="minor", labelsize=tickfontsize)
	#ylim(ymin*.90,ymax*1.10)
	#ylim(1e0,1e2)
	title(titlename)
	subplots_adjust(bottom=0.15)
	savefig(@sprintf("%s_timings_test.png",dataset),bbox_inches="tight")
	savefig(@sprintf("%s_timings_test.pdf",dataset),bbox_inches="tight")

	PyPlot.plt[:show]()
end

showresults(ARGS[1],ARGS[2])

