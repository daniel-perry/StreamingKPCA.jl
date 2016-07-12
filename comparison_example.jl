#! /home/sci/dperry/build/julia/julia/julia

using MAT

using RandomFourierFeatures
using Kernel
using Sketch


X = []
Labels = []

name = "gaussian"
recompute = true

d = 0
n = 0
nt = 0

if name == "gaussian"
 
  n = round(Int64,1e3) 
  d = round(Int64,1e2)
	m = 50 # signal dimension

	eta = 10 # signal to noise ratio

  S = randn(n,d)
  D = eye(d)
  for i=1:m
    D[i,i] = max(1.0 - (i-eta*1.0)/m, 0)
  end
  U,R = qr(randn(d,d))

  X = S * D * U

	if eta > 0
		X += randn(size(X)) / eta
	end
  
	half = round(Int64, ceil(.9*n))
	X[half:end,:] += 1000*randn((n-half+1),d)


	Xtest = X[n/2+1:end,:]
	X = X[1:n/2,:]

	println("X: ", size(X))
	println("Xtest: ", size(Xtest))
 
end

n = size(X,1)
X_centered = X .- (ones(n,1)*mean(X,1))
nt = size(Xtest, 1)
Xtest_centered = Xtest .- (ones(nt,1)*mean(X,1))

function comparison(X, Xtest, K_test, center, sigma)
  n = size(X,1)
  d = size(X,2)

	println("starting comparison")
	FDstyles = ["bo-","ro-","go-","mo-","co-"]
	#ells = [100, 250, 500, 1000, 2000]

	ells = [2,5,10,20,30,50]
	epsilons = 1 ./ ells

	firstval = 10
	lastval = 2000
	incrval = 50

	if n >= 1e4
		firstval = 500
		lastval = 1500
		incrval = 100
	elseif n >= 1e3
		firstval = 500
		lastval = 1500
		incrval = 100
	elseif n >= 1e2
		firstval = 1
		lastval = 200
		incrval = 10
	else
		firstval = 1
		lastval = 100
		incrval = 5
	end

	dims = [firstval:incrval:lastval]
	if(ells[end] < firstval)
		dims = [ells; dims]
	end

	println("dims: ", dims)

	counter = length(dims)

	inputdim = size(X,2)
	errs_fro = Float64[]
	errs_l2 = Float64[]
	timings = Float64[]
	timings_test = Float64[]
	timings_fro = Float64[]
	timings_l2 = Float64[]
	errs_rankk_fro = zeros(length(ells),counter)
	errs_rankk_l2 = zeros(length(ells),counter)
	timings_rankk = zeros(length(ells),counter)
	timings_test_rankk = zeros(length(ells),counter)
	timings_rankk_fro = zeros(length(ells),counter)
	timings_rankk_l2 = zeros(length(ells),counter)
	errs_FD_fro = zeros(length(ells),counter)
	errs_FD_l2 = zeros(length(ells),counter)
	timings_FD = zeros(length(ells),counter)
	timings_test_FD = zeros(length(ells),counter)
	timings_FD_fro = zeros(length(ells),counter)
	timings_FD_l2 = zeros(length(ells),counter)
	errs_nystrom_fro = Float64[]
	errs_nystrom_l2 = Float64[]
	timings_nystrom = Float64[]
	timings_test_nystrom = Float64[]
	timings_nystrom_fro = Float64[]
	timings_nystrom_l2 = Float64[]
	errs_nystromwr_fro = Float64[]
	errs_nystromwr_l2 = Float64[]
	timings_nystromwr = Float64[]
	timings_test_nystromwr = Float64[]
	timings_nystromwr_fro = Float64[]
	timings_nystromwr_l2 = Float64[]
	errs_rff_nystrom_fro = Float64[]
	errs_rff_nystrom_l2 = Float64[]
	timings_rff_nystrom = Float64[]
	timings_test_rff_nystrom = Float64[]
	timings_rff_nystrom_fro = Float64[]
	timings_rff_nystrom_l2 = Float64[]


  rankk_toobig = false
	if n >= 1e2
		rankk_toobig = true
	end
	counter = 0
	for dim=dims #firstval:incrval:lastval
		counter += 1

		samples = 10

		err_fro = zeros(samples)
		err_l2 = zeros(samples)
		timing = zeros(samples)
		timing_test = zeros(samples)
		timing_fro = zeros(samples)
		timing_l2 = zeros(samples)
		err_rankk_fro = zeros(samples,length(ells))
		err_rankk_l2 = zeros(samples,length(ells))
		timing_rankk = zeros(samples,length(ells))
		timing_test_rankk = zeros(samples,length(ells))
		timing_rankk_fro = zeros(samples,length(ells))
		timing_rankk_l2 = zeros(samples,length(ells))
		err_FD_fro = zeros(samples,length(ells))
		err_FD_l2 = zeros(samples,length(ells))
		timing_FD = zeros(samples,length(ells))
		timing_test_FD = zeros(samples,length(ells))
		timing_FD_fro = zeros(samples,length(ells))
		timing_FD_l2 = zeros(samples,length(ells))
		err_nystrom_fro = zeros(samples)
		err_nystrom_l2 = zeros(samples)
		timing_nystrom = zeros(samples)
		timing_test_nystrom = zeros(samples)
		timing_nystrom_fro = zeros(samples)
		timing_nystrom_l2 = zeros(samples)
		err_nystromwr_fro = zeros(samples)
		err_nystromwr_l2 = zeros(samples)
		timing_nystromwr = zeros(samples)
		timing_test_nystromwr = zeros(samples)
		timing_nystromwr_fro = zeros(samples)
		timing_nystromwr_l2 = zeros(samples)
		err_rff_nystrom_fro = zeros(samples)
		err_rff_nystrom_l2 = zeros(samples)
		timing_rff_nystrom = zeros(samples)
		timing_test_rff_nystrom = zeros(samples)
		timing_rff_nystrom_fro = zeros(samples)
		timing_rff_nystrom_l2 = zeros(samples)

		for sample=1:samples
			println( "random fourier features" )
			tic()
			rffs = GenerateFunctionsGaussian(Float64, sigma, inputdim, dim)
			generation_time = toc()

			# time how long it takes to compute K using RFF
			tic()
			mean_compute_time = toc()

			tic()
			Xp = RFFProject( X, rffs )
			if center
				Xp = Xp - ones(size(Xp,1),1) * mean(Xp,1)
			end
			# fair comparison:
			Cp = zeros(size(Xp,2),size(Xp,2))
			fakeweights = ones(size(Xp,1)) 
			for i=1:size(Xp,1)
				Cp += fakeweights[i] * (Xp[i,:]'*Xp[i,:])
			end
			println("compute svd of covariance")
			U,S,V = svd(Cp)

			timing[sample] += toc()
			timing[sample] += generation_time
			timing[sample] += mean_compute_time

			tic()
			Xptest = RFFProject( Xtest, rffs )
			rff_proj_time = toc()

			# full RFF approxmation:
			tic()
			kk = min(ells[end], size(V,2))
			XptestV = Xptest * V
			Kp = XptestV*XptestV'
			timing_test[sample] += toc()
			timing_test[sample] += rff_proj_time

			err_fro[sample] = vecnorm(K_test-Kp) / length(K_test)
			err_l2[sample] = norm(K_test-Kp) / length(K_test[1,:])
	
			for li=1:length(ells)
				l = ells[li]

				println("RFF rank-",l," - projection time")
				tic()
				k = min(dim,l)
				XpV = Xptest * V[:,1:k]
				Kp = XpV * XpV'
				projection_time = toc()
				timing_test_rankk[sample,li] += projection_time

				# time how long it takes to compute the Fro norm
				println( "random fourier features - Fro norm" )
				tic()
				fro_norm_2 = vecnorm(K_test-Kp) / length(K_test)
				err_rankk_fro[sample,li] = fro_norm_2 # frob norm
				timing_rankk_fro[sample,li] += toc() + projection_time + rff_proj_time

				tic()
				l2_norm = norm(K_test-Kp,2) / length(K_test[1,:])
				err_rankk_l2[sample,li] = l2_norm # l2 norm
				timing_rankk_l2[sample,li] += toc() + projection_time + rff_proj_time

				println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
				println("rank-",l," approximation")
				println("sample ", sample, " / ", samples)
				println("rank(K) = ", rank(K_test))
				println("rank(Khat) = ", rank(Kp) )
				println("RFF dim = ", dim, " fro norm = ", fro_norm_2)
				println("RFF dim = ", dim, " l2 norm = ", l2_norm)
				println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			end

  		println( "random fourier features with FD" )
			for li = 1:length(ells)
				l = ells[li]
				p = 1
				if name == "forest"
					p = round(Int64,floor(.5*l))
				end
				println("p = ", p)
				tic()
				println("running FD")
				B = zeros(l,dim)
				B = FD(Xp, l , min(dim,round(Int64, ceil(l/2.0))))
				println("svd of B")
				U,S,V = svd(B)

				timing_FD[sample,li] += toc()
				timing_FD[sample,li] += generation_time
				timing_FD[sample,li] += mean_compute_time

				tic() 
				kk = min(size(V,2),l)
				XpV = Xptest * V[:,1:kk]
				Kp_FD = XpV * XpV'
				timing_test_FD[sample,li] += toc()
				timing_test_FD[sample,li] += rff_proj_time


  			println( "random fourier features with FD - computing fro " )
				tic()
				fro_norm_2 = vecnorm(K_test-Kp_FD) / length(K_test)
				err_FD_fro[sample,li] = fro_norm_2 
				timing_FD_fro[sample,li] += toc()
				tic()
				l2_norm = norm(K_test-Kp_FD,2) / length(K_test[1,:])
				err_FD_l2[sample,li] += l2_norm
				timing_FD_l2[sample,li] += toc()
				println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
				println("sample ", sample, " / ", samples)
				println("RFF dim = ", dim)
				println("rank(K) = ", rank(K_test))
				println("rank(Khat) = ", rank(Kp_FD) )
				println("l = ", l, " fro norm = ", fro_norm_2 )
				println("l = ", l, " l2 norm = ", l2_norm )
				println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			end

			println( "random without replacement sub matrix (Nystrom method)" )
			# reservoir subsample:
			nystrom_num = min(dim, n)  

			UU = zeros(1,1)
			SS = zeros(1,1)
			VV = zeros(1,1)
			X_subset = []
			K_subset = []
			attempts = 0
			while size(UU,1) == 1 && attempts < n
				attempts += 1
				if mod(attempts,1000) == 0
					print(attempts,",")
				end
				# keep sampling until a submatrix that can be factored is found:
				tic()
				subset_ind = collect(1:nystrom_num)
				for j=nystrom_num+1:size(X,1)
					r = round(Int64,ceil(rand()*j))
					if r <= nystrom_num
						# X_subset[r,:] = X[j,:]
						subset_ind[r] = j
					end
				end
				X_subset = X[subset_ind,:]

				K_subset = GaussianKernel(X_subset, sigma)
				try
					UU,SS,VV = svd(K_subset)
				catch(error)
					UU = zeros(1,1)
				end
			end
			K_all = GaussianKernel(Xtest, X_subset, sigma)
			
      epsilon = 1e-10
			if cond(K_subset) > 1e10
        # TSVD:
				SS[ find( SS .<= epsilon ) ] = 0
        ind = find( SS .> epsilon )
        SS[ ind ] = 1 ./ SS[ ind ]
        SS = sqrt(SS)
			else
				SS = sqrt(SS.^(-1))
			end

      epsilon = 1e-10
			if true || cond(K_subset) > 1e10
				UU,SS,VV = svd(K_subset)
        # TSVD:
				k = min(size(UU,2),dim)
        Utkall = UU[:,1:k]' * K_all'
        SS[ find( SS .< epsilon ) ] = 0
        ind = find( SS .> epsilon )
        SS[ ind ] = 1 ./ SS[ ind ]
        SS = sqrt(SS)
        Utkall = diagm(SS) * Utkall
        K_nystrom_wor = Utkall'*Utkall
        # Tikhonov
        #K_subset = UU * diagm(SS .+ epsilon) * VV'
			  #K_nystrom = K_all * (K_subset \ K_all')
      else
			  K_nystrom_wor = K_all * (K_subset \ K_all')
			end
	
			timing_nystrom[sample] += toc()

			println( "random sub matrix (Nystrom method) - computing fro" )
			tic()
			fro_norm_2 = vecnorm(K_test-K_nystrom_wor) / length(K_test)
			err_nystrom_fro[sample] = fro_norm_2
			timing_nystrom_fro[sample] += toc()
			tic()
			l2_norm =  norm(K_test-K_nystrom_wor,2) / length(K_test[1,:])
			err_nystrom_l2[sample] =	l2_norm
			timing_nystrom_l2[sample] += toc()
			println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("sample ", sample, " / ", samples)
			println("RFF dim = ", dim)
			println("rank(K) = ", rank(K_test))
			println("rank(Khat) = ", rank(K_nystrom_wor) )
			println( "nystrom sample size ", size(K_subset,1), ", fro norm = ", fro_norm_2)
			println( "nystrom sample size ", size(K_subset,1), ", l2 norm  = ", l2_norm)
			println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


			println( "random sub matrix with replacement (Nystrom method)" )
			# reservoir subsample:
			nystrom_num = min(dim, n)  

			# fair comparison: iterate over all data points and compute the kernel space norm
			#                  to use importance sampling (reservoir sampling)
			tic()
			X_subset = zeros(nystrom_num,size(X,2))
			norm_est = 0
			for i=nystrom_num+1:size(X,1)
				kern_norm_sq = GaussianKernel(X[i,:],sigma)[1]
				norm_est += kern_norm_sq
				pc = kern_norm_sq / norm_est
				scale = 1 #1/sqrt(nystrom_num*pc)
				for j=1:nystrom_num
					if rand() <= pc
						X_subset[j,:] = scale * X[i,:]
					end
				end
			end

			# do this the same as before to make sure technique above isn't providing biased sample..
			subset_ind = round(Int64, ceil(rand(nystrom_num)*n))
			X_subset = X[subset_ind,:]

			K_subset = GaussianKernel(X_subset, sigma)

			try
				UU,SS,VV = svd(K_subset)
			catch(error)
				UU = zeros(1,1)
			end

			timing_nystromwr[sample] += toc()

			attempts = 0
			while size(UU,1) == 1 && attempts < n # retry the fast way if we picked a bad subset
				attempts += 1
				if mod(attempts,1000) == 0
					print(attempts,",")
				end
				# keep sampling until a submatrix that can be factored is found:
				subset_ind = round(Int64,ceil(rand(nystrom_num)*n))
				X_subset = X[subset_ind,:]

				K_subset = GaussianKernel(X_subset, sigma)
				try
					UU,SS,VV = svd(K_subset)
				catch(error)
					UU = zeros(1,1)
				end
			end

			tic()
			K_all = GaussianKernel(Xtest, X_subset, sigma)

      epsilon = 1e-10
			if cond(K_subset) > 1e10
        # TSVD:
				SS[ find( SS .<= epsilon ) ] = 0
        ind = find( SS .> epsilon )
        SS[ ind ] = 1 ./ SS[ ind ]
        SS = sqrt(SS)
			else
				SS = sqrt(SS.^(-1))
			end

      epsilon = 1e-10
			if true || cond(K_subset) > 1e10
				UU,SS,VV = svd(K_subset)
        # TSVD:
				k = min(size(UU,2),dim)
        Utkall = UU[:,1:k]' * K_all'
        SS[ find( SS .< epsilon ) ] = 0
        ind = find( SS .> epsilon )
        SS[ ind ] = 1 ./ SS[ ind ]
        SS = sqrt(SS)
        Utkall = diagm(SS) * Utkall
        K_nystrom_wr = Utkall'*Utkall
        # Tikhonov
        #K_subset = UU * diagm(SS .+ epsilon) * VV'
			  #K_nystrom = K_all * (K_subset \ K_all')
      else
			  K_nystrom_wr = K_all * (K_subset \ K_all')
			end
	
			timing_test_nystromwr[sample] += toc()

			println( "random sub matrix (Nystrom method) - computing fro" )
			tic()
			fro_norm_2 = vecnorm(K_test-K_nystrom_wr) / length(K_test)
			err_nystromwr_fro[sample] = fro_norm_2
			timing_nystromwr_fro[sample] += toc()
			tic()
			l2_norm =  norm(K_test-K_nystrom_wr,2) / length(K_test[1,:])
			err_nystromwr_l2[sample] =	l2_norm
			timing_nystromwr_l2[sample] += toc()
			println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("sample ", sample, " / ", samples)
			println("rank(K) = ", rank(K_test))
			println("rank(Khat) = ", rank(K_nystrom_wr) )
			println("size(K_subset) = ", size(K_subset))
			println( "nystrom sample size ", size(K_subset,1), ", fro norm = ", fro_norm_2)
			println( "nystrom sample size ", size(K_subset,1), ", l2 norm  = ", l2_norm)
			println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

		end

		push!(errs_fro, median(err_fro))
		push!(errs_l2, median(err_l2))
		push!(timings, mean(timing))
		push!(timings_test, mean(timing_test))

		errs_rankk_fro[:,counter] = median(err_rankk_fro,1)
		errs_rankk_l2[:,counter] = median(err_rankk_l2,1)
		timings_rankk[:,counter] = mean(timing_rankk,1)
		timings_test_rankk[:,counter] = mean(timing_test_rankk,1)
		timings_rankk_fro[:,counter] = mean(timing_rankk_fro,1)
		
		errs_FD_fro[:,counter] = median(err_FD_fro,1)
		errs_FD_l2[:,counter] = median(err_FD_l2,1)
		timings_FD[:,counter] = mean(timing_FD,1)
		timings_test_FD[:,counter] = mean(timing_test_FD,1)
		timings_FD_fro[:,counter] = mean(timing_FD_fro,1)

		push!(errs_nystrom_fro, median(err_nystrom_fro))
		push!(errs_nystrom_l2, median(err_nystrom_l2))
		push!(timings_nystrom, mean(timing_nystrom))
		push!(timings_test_nystrom, mean(timing_test_nystrom))
		push!(timings_nystrom_fro, mean(timing_nystrom_fro))

		push!(errs_nystromwr_fro, median(err_nystromwr_fro))
		push!(errs_nystromwr_l2, median(err_nystromwr_l2))
		push!(timings_nystromwr, mean(timing_nystromwr))
		push!(timings_test_nystromwr, mean(timing_test_nystromwr))
		push!(timings_nystromwr_fro, mean(timing_nystromwr_fro))

		push!(errs_rff_nystrom_fro, median(err_rff_nystrom_fro))
		push!(errs_rff_nystrom_l2, median(err_rff_nystrom_l2))
		push!(timings_rff_nystrom, mean(timing_rff_nystrom))
		push!(timings_test_rff_nystrom, mean(timing_test_rff_nystrom))
		push!(timings_rff_nystrom_fro, mean(timing_rff_nystrom_fro))
		
		#push!(dims, dim)

		println("dim = ", dim, ", FD erro = ", errs_FD_fro[:,counter])
		println("dim = ", dim, ", Nys erro = ", err_nystrom_fro)
		println("dim = ", dim, ", RFF erro = ", err_fro)

		file = matopen(@sprintf("fd_results_%s_n%d_nt%d_d%d_s%f.mat",name,n,nt,d,sigma),"w")
		write(file, "dims", dims)
		write(file, "ells", ells)
		write(file, "epsilons", epsilons)
		write(file, "errs_fro", errs_fro)
		write(file, "errs_l2", errs_l2)
		write(file, "timings", timings)
		write(file, "timings_test", timings_test)
		write(file, "timings_fro", timings_fro)
		write(file, "timings_l2", timings_l2)
		write(file, "errs_rankk_fro", errs_rankk_fro)
		write(file, "errs_rankk_l2", errs_rankk_l2)
		write(file, "timings_rankk", timings_rankk)
		write(file, "timings_test_rankk", timings_test_rankk)
		write(file, "timings_rankk_fro", timings_rankk_fro)
		write(file, "timings_rankk_l2", timings_rankk_l2)
		write(file, "errs_FD_fro", errs_FD_fro)
		write(file, "errs_FD_l2", errs_FD_l2)
		write(file, "timings_FD", timings_FD)
		write(file, "timings_test_FD", timings_test_FD)
		write(file, "timings_FD_fro", timings_FD_fro)
		write(file, "timings_FD_l2", timings_FD_l2)
		write(file, "errs_nystrom_fro", errs_nystrom_fro)
		write(file, "errs_nystrom_l2", errs_nystrom_l2)
		write(file, "timings_nystrom", timings_nystrom)
		write(file, "timings_test_nystrom", timings_test_nystrom)
		write(file, "timings_nystrom_fro", timings_nystrom_fro)
		write(file, "timings_nystrom_l2", timings_nystrom_l2)
		write(file, "errs_nystromwr_fro", errs_nystromwr_fro)
		write(file, "errs_nystromwr_l2", errs_nystromwr_l2)
		write(file, "timings_nystromwr", timings_nystromwr)
		write(file, "timings_test_nystromwr", timings_test_nystromwr)
		write(file, "timings_nystromwr_fro", timings_nystromwr_fro)
		write(file, "timings_nystromwr_l2", timings_nystromwr_l2)
		write(file, "errs_rff_nystrom_fro", errs_rff_nystrom_fro)
		write(file, "errs_rff_nystrom_l2", errs_rff_nystrom_l2)
		write(file, "timings_rff_nystrom", timings_rff_nystrom)
		write(file, "timings_test_rff_nystrom", timings_test_rff_nystrom)
		write(file, "timings_rff_nystrom_fro", timings_rff_nystrom_fro)
		write(file, "timings_rff_nystrom_l2", timings_rff_nystrom_l2)
		close(file)
	end

	file = matopen(@sprintf("fd_results_%s_n%d_nt%d_d%d_s%f.mat",name,n,nt,d,sigma),"w")
	write(file, "timing_orig", timing_orig)
	write(file, "dims", dims)
	write(file, "ells", ells)
	write(file, "epsilons", epsilons)
	write(file, "errs_fro", errs_fro)
	write(file, "errs_l2", errs_l2)
	write(file, "timings", timings)
	write(file, "timings_test", timings_test)
	write(file, "timings_fro", timings_fro)
	write(file, "timings_l2", timings_l2)
	write(file, "errs_rankk_fro", errs_rankk_fro)
	write(file, "errs_rankk_l2", errs_rankk_l2)
	write(file, "timings_rankk", timings_rankk)
	write(file, "timings_test_rankk", timings_test_rankk)
	write(file, "timings_rankk_fro", timings_rankk_fro)
	write(file, "timings_rankk_l2", timings_rankk_l2)
	write(file, "errs_FD_fro", errs_FD_fro)
	write(file, "errs_FD_l2", errs_FD_l2)
	write(file, "timings_FD", timings_FD)
	write(file, "timings_test_FD", timings_test_FD)
	write(file, "timings_FD_fro", timings_FD_fro)
	write(file, "timings_FD_l2", timings_FD_l2)
	write(file, "errs_nystrom_fro", errs_nystrom_fro)
	write(file, "errs_nystrom_l2", errs_nystrom_l2)
	write(file, "timings_nystrom", timings_nystrom)
	write(file, "timings_test_nystrom", timings_test_nystrom)
	write(file, "timings_nystrom_fro", timings_nystrom_fro)
	write(file, "timings_nystrom_l2", timings_nystrom_l2)
	write(file, "errs_nystromwr_fro", errs_nystromwr_fro)
	write(file, "errs_nystromwr_l2", errs_nystromwr_l2)
	write(file, "timings_nystromwr", timings_nystromwr)
	write(file, "timings_test_nystromwr", timings_test_nystromwr)
	write(file, "timings_nystromwr_fro", timings_nystromwr_fro)
	write(file, "timings_nystromwr_l2", timings_nystromwr_l2)
	write(file, "errs_rff_nystrom_fro", errs_rff_nystrom_fro)
	write(file, "errs_rff_nystrom_l2", errs_rff_nystrom_l2)
	write(file, "timings_rff_nystrom", timings_rff_nystrom)
	write(file, "timings_test_rff_nystrom", timings_test_rff_nystrom)
	write(file, "timings_rff_nystrom_fro", timings_rff_nystrom_fro)
	write(file, "timings_rff_nystrom_l2", timings_rff_nystrom_l2)
	close(file)
	
	#PyPlot.plt.show()

end

centerit = false
K = []
Ktest = []
timing_orig = 0
if recompute
	tic()

	# compute sigma from data:
	sigma = 0
	println("computing sigma parameter...")
	if n > 1e3
		ss = round(Int64,1e3)
		sampled = X[1:ss,:]
		for i=ss+1:n
			rr = round(Int64,ceil(rand()*n))
			if rr <= ss
				sampled[rr,:] = X[i,:]
			end
		end
		for i=1:ss
			sigma += sum(sqrt(sum((sampled - ones(ss,1)*sampled[i,:]).^2,2)))
		end
		sigma /= ss^2
	else
		for i=1:n
			sigma += sum(sqrt(sum((X - ones(n,1)*X[i,:]).^2,2)))
		end
		sigma /= n^2
	end

	println("computed sigma: ", sigma)

	println("computing kernel matrix...")
	K_test = GaussianKernel(Xtest_centered,sigma)
	nn = size(K,1)
	if centerit
		println("centering kernel matrix...")
		K_test = K_test - (1/nn)*ones(size(K))*K - K*(1/nn)*ones(size(K)) + (1/nn)*ones(size(K))*K*(1/nn)*ones(size(K))
	end
	timing_orig = toc()

	# write out K and timing:
	file = matopen(@sprintf("fd_results_K_%s_%d.mat",name,n),"w")
	write(file,"X",X)
	write(file,"X_centered",X_centered)
	write(file,"Xtest",Xtest)
	write(file,"Xtest_centered",Xtest_centered)
	write(file,"K_test",K_test)
	write(file,"sigma",sigma)
	write(file,"timing_orig",timing_orig)
	close(file)
else
	file = matopen(@sprintf("fd_results_K_%s_%d.mat",name,n))
	K_test = read(file,"K_test")
	X = read(file,"X")
	X_centered = read(file,"X_centered")
	Xtest = read(file,"Xtest")
	Xtest_centered = read(file,"Xtest_centered")
	sigma = read(file,"sigma")
	timing_orig = read(file,"timing_orig")
	println("timing to compute K: ", timing_orig)
	close(file)
end

comparison(X_centered, Xtest_centered, K_test, centerit, sigma)
