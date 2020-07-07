# Joulosky Airfoil Camberline Module
# Developed By Yi Tsung Lee (Joe)
# Version  V 1.1 (0622 2020)

#====== incase no SymEngine Package installed =========
Pkg.add("SymEngine")
=#
using Plots
gr()
using LinearAlgebra
# using SymEngine
# using PotentialFlow

function Newton_method(f::Function,f0::Function,x0::Float64,
  tol::Float64 = 1e-5,maxiter::Integer=200, eps::Float64=1e-10,m::Float64 = 0.7)

     for i = 1:maxiter
       #print("iter =",i,"\n")
       f0_value = f0(x0)
       if abs(f0_value) < eps
           @warn("first derivative is zero! \n")
          return x0
       end
       f_value = f(x0)
       x_n1 = x0 - m*f_value/f0_value
       if abs((x_n1-x0)/x0) < tol
          #print("sol =",x_n1,"\n")
          return x_n1
       end
       x0 = x_n1
     end
    print("sol =",x_n1,"\n")
    @error("max iteration reached \n")
end

mutable struct circle_data
    k ::Float64
    # singularity point of circle
    xi_s  ::Float64
    eta_s ::Float64
    # center point of circle
    xi_c  ::Float64
    eta_c ::Float64
    # initial circle radius
    a ::Float64

    # Circle  coordinate space
    theta ::Array{Float64}
    xi ::Array{Float64}
    eta::Array{Float64}
    theta_TE::Float64
    theta_LE::Float64
    xi_TE::Float64
    eta_TE::Float64
    xi_LE::Float64
    eta_LE::Float64

    function circle_data(k,xi_s,eta_s,xi_c,eta_c)
       k = k
       xi_s  = xi_s
       eta_s = eta_s
       xi_c  = xi_c
       eta_c = eta_c
       new(k,xi_s,eta_s,xi_c,eta_c)
    end
end

function circle_generation(data::circle_data,n_intp::Int64)
    k = data.k
    xi_s  = data.xi_s
    eta_s = data.eta_s
    xi_c  = data.xi_c
    eta_c = data.eta_c

    data.a = sqrt((xi_s-xi_c)^2+(eta_s-eta_c)^2)
    data.theta_TE = mod(2*pi+atan((eta_s-eta_c)/(xi_s-xi_c)),2*pi)
    data.theta_LE = mod(pi + -1*atan((eta_s-eta_c)/(xi_s-xi_c)),2*pi)

    zeta = xi_c + eta_c*im + k*data.a*exp(data.theta_TE*im)
    data.xi_TE = real(zeta)
    data.eta_TE = imag(zeta)
    zeta = xi_c + eta_c*im + k*data.a*exp(data.theta_LE*im)
    data.xi_LE = real(zeta)
    data.eta_LE = imag(zeta)

    # circle coordinates
    data.theta = zeros(n_intp*2+2)
    data.theta[1] = data.theta_TE
    data.theta[2+n_intp] = data.theta_LE
    lo_range = (data.theta_TE-data.theta_LE)
    up_range = 2*pi-lo_range
    for i = 1:n_intp
     data.theta[i+1] = mod(i*up_range/(n_intp+1)+data.theta_TE,2*pi)
     data.theta[i+2+n_intp] = mod(i*lo_range/(n_intp+1)+data.theta_LE,2*pi)
    end

    data.xi = zeros(n_intp*2+2)
    data.eta = zeros(n_intp*2+2)
    for i = 1:n_intp*2+2
     zeta = xi_c + eta_c*im + k*data.a*exp(data.theta[i]*im)
     data.xi[i] = real(zeta)
     data.eta[i] = imag(zeta)
    end

   return data
end

mutable struct airfoil_data

      n_intp::Int64
      # airfoil coordinate
      x::Array{Float64}
      y::Array{Float64}
      t::Array{Float64}

      # camber coordinate
      c_x::Array{Float64}
      c_y::Array{Float64}

      #

      function airfoil_data(n_intp::Int64)
      n_inttp = n_intp
      new(n_intp)
     end
end

function joukowsky(cd::circle_data,ad::airfoil_data)
   z = cd.xi + cd.eta.*im + (cd.xi_s+cd.eta_s.*im).^2 ./ (cd.xi + cd.eta.*im)
   ad.x = real(z)
   ad.y = imag(z)

   return ad
end

num_intp = 100

data = circle_data(1.2,1.0,0.0,-0.2,0.2)
data = circle_generation(data,num_intp)
plot(data.xi,data.eta,seriestype = :scatter,aspect_ratio =:equal,size = (1200, 800),dpi = 300)
plot!([data.xi[1],data.xi_c],[data.eta[1],data.eta_c],markershape = :circle)
plot!([data.xi[num_intp+2],data.xi_c],[data.eta[num_intp+2],data.eta_c],markershape = :circle)
savefig("circle.png")
#=
function num_camber(cd::circle_data,ad::airfoil_data,iter_limit::Int64)

   # ===== Start of initializationin the method =========
   # initialize camber point data set in ad
   ad.c_x = zeros(Float64,ad.n_intp+2)
   ad.c_y = zeros(Float64,ad.n_intp+2)

   # Assign TE/LE point
   ad.c_x[1] = ad.x[1]
   ad.c_y[1] = ad.y[1]
   ad.c_x[ad.n_intp+2] = ad.x[ad.n_intp+2]
   ad.c_y[ad.n_intp+2] = ad.y[ad.n_intp+2]

   # Assign interpolation mid camberpoint
    for i = 1:ad.n_intp
     ad.c_x[1+i] = (ad.x[1+i] + ad.x[end+1-i])/2
     ad.c_y[1+i] = (ad.y[1+i] + ad.y[end+1-i])/2
    end

    # starting initializing the slope of numerical camber line

    slp_1 = zeros(Float64,ad.n_intp+1) # camber sectional slope
    slp_2 = zeros(Float64,ad.n_intp)   # camber sectional average slope
    nom = zeros(Float64,ad.n_intp)     # camber sectional average normal
    c_x_old = zeros(Float64,ad.n_intp + 2)
    c_y_old = zeros(Float64,ad.n_intp + 2)

    # initialize the theta matrix
    numtheta_up = zeros(Float64,ad.n_intp)
    numtheta_lo = zeros(Float64,ad.n_intp)
    numtheta_TE = cd.theta_TE
    numtheta_LE = cd.theta_LE
    #print("initial num_theta_LE = ",numtheta_LE,"\n")
    #print("initial num_theta_TE = ",numtheta_TE,"\n")
    for i = 1:ad.n_intp
      numtheta_up[i] = cd.theta[i+1]  # store from camber point 1->n_intp
      numtheta_lo[i] = cd.theta[end+1-i] # store from camber point 1->n_intp
    end

    # predefine function of coordinate needed in the method
    xi(theta)  = cd.xi_c +cd.k*cd.a*cos(theta)
    eta(theta) = cd.eta_c +cd.k*cd.a*sin(theta)

    r(theta)   = sqrt(xi(theta)^2 + eta(theta)^2)
    dr(theta)  = (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))

                 #(2*sqrt(cd.xi_c^2 + cd.eta_c^2 + (cd.k*cd.a)^2 +
                 #2*cd.xi_c*cd.k*cd.a*cos(theta) + 2*cd.eta_c*cd.k*cd.a*sin(theta) ))
    dr2(theta) = (-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2/r(theta)^3

    dr3(theta) = (2*cd.xi_c*cd.k*cd.a*sin(theta) - 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^3-
                 1/2*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*cos(theta)-2*cd.eta_c*cd.k*cd.a*sin(theta))/r(theta)^3+
                 3/8*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^5


    x(theta)  = xi(theta)*(1 + 1/r(theta)^2)
    dx(theta) = -1*cd.k*cd.a*sin(theta)*(1 + 1/r(theta)^2) +
                xi(theta)*(-2/r(theta)^3)*dr(theta)
    dx2(theta) = xi(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*sin(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*cos(theta)*(1+1/r(theta)^2)
    dx3(theta) = xi(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*sin(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 + 1 + 1/r(theta)^2)+
                 cd.k*cd.a*cos(theta)*(6*dr(theta)/r(theta)^3)


    y(theta)  = eta(theta)*(1 - 1/r(theta)^2)
    dy(theta) = cd.k*cd.a*cos(theta)*(1 - 1/r(theta)^2) +
                eta(theta)*(2/r(theta)^3)*dr(theta)
    dy2(theta) = -1*eta(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*cos(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*sin(theta)*(1-1/r(theta)^2)
    dy3(theta) = -1*eta(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*cos(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 - 1 + 1/r(theta)^2)-
                 cd.k*cd.a*sin(theta)*(6*dr(theta)/r(theta)^3)

     # curvature module
    kappa(theta) = (dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)
    dkappa(theta) = (dx(theta)*dy3(theta)-dy(theta)*dx3(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)-
                    3/2*(dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(5/2)*(2*dx(theta)*dx2(theta)+2*dy(theta)*dy2(theta))

    total_point = 2*ad.n_intp+2
    kpa = zeros(Float64,total_point)
    phi = zeros(Float64,total_point)
    phi_edge = zeros(Float64,2)
    edge_index = 1
    for i = 1:total_point
        phi[i] = 2*pi/total_point*(i-1)
        kpa[i] = kappa(phi[i])
    end

    # find the local max point
    for i = 1:total_point
        L = i-1
        R = i+1
      if i == 1
        L = total_point
      elseif i == total_point
        R = 1
      end

      if kpa[i]-kpa[L]>=0 && kpa[i]-kpa[R]>=0
        phi_edge[edge_index] = phi[i]
        edge_index = edge_index + 1
        if edge_index == 3
            break
         end
       end
    end

    #plot(phi,kpa,xlabel ="theta",ylabel = "/kappa",legend = false)
    #savefig("kappa.png")
    print("extension LE =",numtheta_LE,"\n")
    print("kappa LE =",phi_edge[1],"\n")
    print("extension TE =",numtheta_TE,"\n")
    print("kappa TE =",phi_edge[2],"\n")
    # =========== end of kappa module ==================

    # =========== use kappa TE/LE to redefine initial theta point
    numtheta_TE = phi_edge[2]
    numtheta_LE = phi_edge[1]
    # use updated LE/TE position to redefine airfoil surface theta
    lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
    up_range = 2*pi-lo_range

    for i = 1:ad.n_intp
     numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
     numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
     # after renewing the TE/LE , recalculate the initial camberline
     ad.c_x[i+1] = (x(numtheta_up[i]) + x(numtheta_lo[i]))/2
     ad.c_y[i+1] = (y(numtheta_up[i]) + y(numtheta_lo[i]))/2
     x_up = x(numtheta_up[i])
     y_up = y(numtheta_up[i])
     x_lo = x(numtheta_lo[i])
     y_lo = y(numtheta_lo[i])
     plot!([x_lo,x_up],[y_lo,y_up],legend = false)
    end

    ad.c_x[1] = x(numtheta_TE)
    ad.c_y[1] = y(numtheta_TE)
    ad.c_x[ad.n_intp+2] = x(numtheta_LE)
    ad.c_y[ad.n_intp+2] = y(numtheta_LE)
    # =========== end of use kappa TE/LE to redefine initial theta point

    # ===== end of initializationin the method =========

    # ===== Start of the iteration loop ================
    for iter = 1:iter_limit

       # step 1. extract all section slope and normal value from current camber
       # save the current camber point into x_c_old, y_c_old
       for i = 1: ad.n_intp + 2
         c_x_old[i] = ad.c_x[i]
         c_y_old[i] = ad.c_y[i]
       end

       # calculate the numerical camber lope
       for i = 1:ad.n_intp +1
         slp_1[i] = (ad.c_y[i+1]- ad.c_y[i])/(ad.c_x[i+1]- ad.c_x[i])
       end
       #  6/22 update : using weighted slope calculation (weighting from camber section length)
       for i = 1:ad.n_intp
         l1 = sqrt((ad.c_x[i+1]-ad.c_x[i])^2+(ad.c_y[i+1]-ad.c_y[i])^2)
         l2 = sqrt((ad.c_x[i+2]-ad.c_x[i+1])^2+(ad.c_y[i+2]-ad.c_y[i+1])^2)
         slp_2[i] = (slp_1[i]*l1 + slp_1[i+1]*l2)/(l1+l2)
         nom[i] = -1/slp_2[i]
       end
       # 6/17 update , unhook the connection of the LE/TE slope for the first
       # and last camber point(TE/LE point)

       # unhook doesn't mean that you should only use the second slope .Figure out another way
       # 6/22 idea : Once renew the LE/TE ,renew the initial camberline, since the initial camberline
       # will determine even the second slope is usable as approximation
       #=
       #  6/22 update: after each iteration, update of the TE and LE might be needed
       # ================  TE/LE shifting module =========================
       if slp_1[1]*slp_1[2] < 0   # TE need to redefine
           nom[1] = -1/slp_1[2]
           local A = slp_1[2]
           local B = ad.c_y[2] - A*ad.c_x[2]

           f1(theta) = y(theta) - A*x(theta)-B
           df1(theta) = dy(theta) - A*dx(theta)
           numtheta_TE = mod(Newton_method(f1,df1,numtheta_TE),2*pi)
           ad.c_x[1] = x(numtheta_TE)
           ad.c_y[1] = y(numtheta_TE)

           print("iter=$iter TE Modification finished\n")
           print("new num_theta_TE = ",numtheta_TE,"\n")
       end
       if slp_1[ad.n_intp + 1]*slp_1[ad.n_intp] < 0 # lE need to redefine
           nom[ad.n_intp] = -1/slp_1[ad.n_intp]
           local A = slp_1[ad.n_intp]
           local B = ad.c_y[ad.n_intp + 1] - A*ad.c_x[ad.n_intp + 1]

           f2(theta) = y(theta) - A*x(theta)-B
           df2(theta) = dy(theta) - A*dx(theta)

           numtheta_LE = mod(Newton_method(f2,df2,numtheta_LE),2*pi)
           ad.c_x[ad.n_intp+2] = x(numtheta_LE)
           ad.c_y[ad.n_intp+2] = y(numtheta_LE)
           print("iter=$iter LE Modification finished\n")
           print("new num_theta_LE = ",numtheta_LE,"\n")
       end

       # use updated LE/TE position to redefine airfoil surface theta
       lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
       up_range = 2*pi-lo_range
       for i = 1:ad.n_intp
        numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
        numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
        # after renewing the TE/LE , recalculate the initial camberline
        ad.c_x[i+1] = (x(numtheta_up[i]) + x(numtheta_lo[i]))/2
        ad.c_y[i+1] = (y(numtheta_up[i]) + y(numtheta_lo[i]))/2
       end

       # ================  end of TE/LE shifting module =========================
       =#
       # end of step 1.

      # step 2 , use sectional normal and the camber point to create
      # intersection line y = Ax+B

      # loop through all camber point (interpolation camber point)

    for i = 1 : ad.n_intp

         local A = nom[i]
         local B = ad.c_y[i+1] - A*ad.c_x[i+1]

         f(theta) = y(theta) - A*x(theta)-B
         df(theta) = dy(theta) - A*dx(theta)

      # set up the initial guess using the theta angle correponding to i position
      # on both up and lower surface

      iden = 1 # flag to identify the calculationcondition
      #=================================================
       iden = 1: upper theta in lower surface
       iden = 2: lower theta in upper surface
       iden = 0: iteration result looks reasonable
      ==================================================#
      theta_up0 = numtheta_up[i] # initial guess value
      theta_lo0 = numtheta_lo[i] # initial guess value
      theta_up = 0 # upper intersection point theta
      theta_lo = 0 # lower intersection point theta

      shift = 0 # index for adjusting initial guess

     # examine and initial guess shifting module
     while (iden == 1||iden ==2)
       print(" iter = $iter at camberpoint $i , initial guess shiting =",shift,"\n")
       shift = shift + 1
       # step 3 solve the equatio
       theta_up = mod(Newton_method(f,df,theta_up0),2*pi)
       theta_lo = mod(Newton_method(f,df,theta_lo0),2*pi)

      # send warning if the iteration method fail to capture 2 different theta
      # on upper an dlower surface
      if theta_up == theta_lo # this condition is nearly impossible
          @warn "fail to capture two different angle at $(i+1) th camber point"
          iden = 0

      elseif theta_up >= numtheta_TE &&
             theta_up <= numtheta_LE
          @warn " $(i) th camber point has upper theta on the lower surface !!"
          iden = 1
          # near LE, use previous theta
          if i-shift > 0
              theta_up0 = numtheta_up[i-shift]
          # near TE, use next theta
          elseif i-shift <0
              theta_up0 = numtheta_up[i+shift]
          end
      elseif theta_lo >= numtheta_TE ||
             theta_lo <= numtheta_LE
          @warn " $(i) th camber point has lower theta on the upper surface !!"
          iden = 2
          # near LE, use previous theta
          if i-shift > 0
              theta_lo0 = numtheta_lo[i-shift]
          # near TE, use next theta
          elseif i-shift <0
              theta_lo0 = numtheta_lo[i+shift]
          end

      else #no problem
          print("iteration success, current camber point  theta are in the range \n")
          iden = 0
      end
     end
     # end of initial guess shifting module

      x_up = x(theta_up)
      y_up = y(theta_up)
      x_lo = x(theta_lo)
      y_lo = y(theta_lo)

      # if it is the last iteration, mark the intersection line
      if  iter == iter_limit
         plot!([x_lo,x_up],[y_lo,y_up],legend = false)
      end

      # update the new camber line interpolation data

        examinexy[i,1:4] = [x_up,y_up,x_lo,y_lo]
        ad.c_x[i+1] = (x_up + x_lo)/2
        ad.c_y[i+1] = (y_up + y_lo)/2
      end
      # end of 1 whoel iteration to create a line camberline
      # print(slp_1[ad.n_intp + 1])

      # juedgement for the convergence :
      # average change of distance of the iteration on single camber point
      # to the chord length smaller
      chord  = sqrt((ad.c_x[1]-ad.c_x[end])^2 + (ad.c_y[1]-ad.c_y[end])^2)
      ds = 0
       for i = 1: ad.n_intp+2
          ds = ds + sqrt(((ad.c_x[i]-c_x_old[i])^2 + (ad.c_y[i]-c_y_old[i])^2)/ad.n_intp)
       end
      if ds/chord <= 1e-5
         print("camberline iter finished \n")
         break # break itertaion loop
      end

    end
   return ad
end
function num_camber2(cd::circle_data,ad::airfoil_data,iter_limit::Int64)

   # ===== Start of initializationin the method =========
   # initialize camber point data set in ad
   ad.c_x = zeros(Float64,ad.n_intp+2)
   ad.c_y = zeros(Float64,ad.n_intp+2)

    # starting initializing the slope of numerical camber line

    slp_1 = zeros(Float64,ad.n_intp+1) # camber sectional slope
    slp_2 = zeros(Float64,ad.n_intp)   # camber sectional average slope
    nom = zeros(Float64,ad.n_intp)     # camber sectional average normal
    c_x_old = zeros(Float64,ad.n_intp + 2)
    c_y_old = zeros(Float64,ad.n_intp + 2)

    # initialize the theta matrix
    numtheta_up = zeros(Float64,ad.n_intp)
    numtheta_lo = zeros(Float64,ad.n_intp)
    numtheta_TE = 0.0
    numtheta_LE = 0.0

    # predefine function of coordinate needed in the method
    xi(theta)  = cd.xi_c +cd.k*cd.a*cos(theta)
    eta(theta) = cd.eta_c +cd.k*cd.a*sin(theta)

    r(theta)   = sqrt(xi(theta)^2 + eta(theta)^2)
    dr(theta)  = (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))

    dr2(theta) = (-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2/r(theta)^3

    dr3(theta) = (2*cd.xi_c*cd.k*cd.a*sin(theta) - 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^3-
                 1/2*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*cos(theta)-2*cd.eta_c*cd.k*cd.a*sin(theta))/r(theta)^3+
                 3/8*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^5


    x(theta)  = xi(theta)*(1 + 1/r(theta)^2)
    dx(theta) = -1*cd.k*cd.a*sin(theta)*(1 + 1/r(theta)^2) +
                xi(theta)*(-2/r(theta)^3)*dr(theta)
    dx2(theta) = xi(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*sin(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*cos(theta)*(1+1/r(theta)^2)
    dx3(theta) = xi(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*sin(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 + 1 + 1/r(theta)^2)+
                 cd.k*cd.a*cos(theta)*(6*dr(theta)/r(theta)^3)


    y(theta)  = eta(theta)*(1 - 1/r(theta)^2)
    dy(theta) = cd.k*cd.a*cos(theta)*(1 - 1/r(theta)^2) +
                eta(theta)*(2/r(theta)^3)*dr(theta)
    dy2(theta) = -1*eta(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*cos(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*sin(theta)*(1-1/r(theta)^2)
    dy3(theta) = -1*eta(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*cos(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 - 1 + 1/r(theta)^2)-
                 cd.k*cd.a*sin(theta)*(6*dr(theta)/r(theta)^3)

    # =================== curvature(kappa) module ==========================
    kappa(theta) = (dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)
    dkappa(theta) = (dx(theta)*dy3(theta)-dy(theta)*dx3(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)-
                    3/2*(dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(5/2)*(2*dx(theta)*dx2(theta)+2*dy(theta)*dy2(theta))

    total_point = 2*ad.n_intp+2
    kpa = zeros(Float64,total_point)
    phi = zeros(Float64,total_point)
    phi_edge = zeros(Float64,2)
    edge_index = 1
    for i = 1:total_point
        phi[i] = 2*pi/total_point*(i-1)
        kpa[i] = kappa(phi[i])
    end

    # find the local max point
    for i = 1:total_point
        L = i-1
        R = i+1
      if i == 1
        L = total_point
      elseif i == total_point
        R = 1
      end

      if kpa[i]-kpa[L]>=0 && kpa[i]-kpa[R]>=0
        phi_edge[edge_index] = phi[i]
        edge_index = edge_index + 1
        if edge_index == 3
            break
         end
       end
    end

    #plot(phi,kpa,xlabel ="theta",ylabel = "/kappa",legend = false)
    #savefig("kappa.png")
    print("kappa LE theta =",phi_edge[1],"\n")
    print("kappa TE theta =",phi_edge[2],"\n")
    # =========== end of kappa module ==================
    # redefine data at LE/TE
    numtheta_TE = phi_edge[2]
    numtheta_LE = phi_edge[1]
    ad.c_x[1] = x(numtheta_TE)
    ad.c_y[1] = y(numtheta_TE)
    ad.c_x[ad.n_intp+2] = x(numtheta_LE)
    ad.c_y[ad.n_intp+2] = y(numtheta_LE)
    ad.x[1] = ad.c_x[1]
    ad.y[1] = ad.c_y[1]
    ad.x[ad.n_intp+2] = ad.c_x[ad.n_intp+2]
    ad.y[ad.n_intp+2] = ad.c_y[ad.n_intp+2]
    slp_1[1] = dy(numtheta_TE)/dx(numtheta_TE)
    slp_1[ad.n_intp] = dy(numtheta_LE)/dx(numtheta_LE)


    # use updated LE/TE position to redefine airfoil surface theta
    # these value will be used as initial guess for NEWTON Method
    lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
    up_range = 2*pi-lo_range
    for i = 1:ad.n_intp
     numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
     numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
     ad.c_x[i+1] = (x(numtheta_up[i]) + x(numtheta_lo[i]))/2
     ad.c_y[i+1] = (y(numtheta_up[i]) + y(numtheta_lo[i]))/2
     ad.x[i+1] = x(numtheta_up[i])
     ad.y[i+1] = y(numtheta_up[i])
     ad.x[end+1-i] = x(numtheta_lo[i])
     ad.y[end+1-i] = y(numtheta_lo[i])
    end
    # ===== end of initializationin the method =========

    #=
    # generate the initial camber using kappa
    for i = 1:ad.n_intp
        # find the pairing for upper point
        k_theta_up = numtheta_up[i]
        k = kappa(k_theta_up)
        f(theta) = kappa(theta)-k
        df(theta) = dkappa(theta)
        k_theta_lo = mod(Newton_method(f,df,numtheta_lo[i]),2*pi)

        x_up = x(k_theta_up)
        y_up = y(k_theta_up)
        x_lo = x(k_theta_lo)
        y_lo = y(k_theta_lo)
        ad.c_x[i+1] = (x_up + x_lo)/2
        ad.c_y[i+1] = (y_up + y_lo)/2

        plot!([x_lo,x_up],[y_lo,y_up],legend = false)

    end
    =#

    # ===== Start of the iteration loop ================
    for iter = 1:iter_limit

       # step 1. extract all section slope and normal value from current camber
       # save the current camber point into x_c_old, y_c_old
       for i = 1: ad.n_intp + 2
         c_x_old[i] = ad.c_x[i]
         c_y_old[i] = ad.c_y[i]
       end

       # calculate the numerical camber lope
       for i = 1:ad.n_intp +1
         slp_1[i] = (ad.c_y[i+1]- ad.c_y[i])/(ad.c_x[i+1]- ad.c_x[i])
       end

       #6/28 update : Force the edge slope to be
       slp_1[1] = -1*dx(numtheta_TE)/dy(numtheta_TE)
       slp_1[end] = -1*dx(numtheta_LE)/dy(numtheta_LE)

       #  6/22 update : using weighted slope calculation (weighting from camber section length)
       for i = 1:ad.n_intp
         l1 = sqrt((ad.c_x[i+1]-ad.c_x[i])^2+(ad.c_y[i+1]-ad.c_y[i])^2)
         l2 = sqrt((ad.c_x[i+2]-ad.c_x[i+1])^2+(ad.c_y[i+2]-ad.c_y[i+1])^2)
         slp_2[i] = (slp_1[i]*l1 + slp_1[i+1]*l2)/(l1+l2)
         nom[i] = -1/slp_2[i]
       end
       # for reverse slope , jump to use average
       for i = 1:ad.n_intp
         if slp_1[i]*slp_1[i+1] <=0
           slp_2[i] = (ad.c_y[i+2]- ad.c_y[i])/(ad.c_x[i+2]- ad.c_x[i])
           nom[i] = -1/slp_2[i]
         end
       end

       # 6/28 update : fix the edge slope as  perpendicular to the edge slope
       nom[1] = dy(numtheta_TE)/dx(numtheta_TE)
       nom[ad.n_intp] = dy(numtheta_LE)/dx(numtheta_LE)

       # 6/17 update , unhook the connection of the LE/TE slope for the first
       # and last camber point(TE/LE point)

       # unhook doesn't mean that you should only use the second slope .Figure out another way
       # 6/22 idea : Once renew the LE/TE ,renew the initial camberline, since the initial camberline
       # will determine even the second slope is usable as approximation


       # end of step 1.

      # step 2 , use sectional normal and the camber point to create
      # intersection line y = Ax+B

      # loop through all camber point (interpolation camber point)

    for i = 1 : ad.n_intp

         local A = nom[i]
         local B = ad.c_y[i+1] - A*ad.c_x[i+1]

         f(theta) = y(theta) - A*x(theta)-B
         df(theta) = dy(theta) - A*dx(theta)

      # set up the initial guess using the theta angle correponding to i position
      # on both up and lower surface

      iden = 1 # flag to identify the calculationcondition
      #=================================================
       iden = 1: upper theta in lower surface
       iden = 2: lower theta in upper surface
       iden = 0: iteration result looks reasonable
      ==================================================#
      theta_up0 = numtheta_up[i] # initial guess value
      theta_lo0 = numtheta_lo[i] # initial guess value
      theta_up = 0 # upper intersection point theta
      theta_lo = 0 # lower intersection point theta

      shift = 0 # index for adjusting initial guess

     # examine and initial guess shifting module
     while (iden == 1||iden ==2)
       print(" iter = $iter at camberpoint $i , initial guess shiting =",shift,"\n")
       shift = shift + 1
       # step 3 solve the equatio
       theta_up = mod(Newton_method(f,df,theta_up0),2*pi)
       theta_lo = mod(Newton_method(f,df,theta_lo0),2*pi)

      # send warning if the iteration method fail to capture 2 different theta
      # on upper an dlower surface
      if theta_up == theta_lo # this condition is nearly impossible
          @warn "fail to capture two different angle at $(i+1) th camber point"
          iden = 0

      elseif theta_up >= numtheta_TE &&
             theta_up <= numtheta_LE
          @warn " $(i) th camber point has upper theta on the lower surface !!"
          iden = 1
          # near LE, use previous theta
          if i-shift > 0
              theta_up0 = numtheta_up[i-shift]
          # near TE, use next theta
          elseif i-shift <0
              theta_up0 = numtheta_up[i+shift]
          end
      elseif theta_lo >= numtheta_TE ||
             theta_lo <= numtheta_LE
          @warn " $(i) th camber point has lower theta on the upper surface !!"
          iden = 2
          # near LE, use previous theta
          if i-shift > 0
              theta_lo0 = numtheta_lo[i-shift]
          # near TE, use next theta
          elseif i-shift <0
              theta_lo0 = numtheta_lo[i+shift]
          end

      else #no problem
          print("iteration success, current camber point  theta are in the range \n")
          iden = 0
      end
     end
     # end of initial guess shifting module

      x_up = x(theta_up)
      y_up = y(theta_up)
      x_lo = x(theta_lo)
      y_lo = y(theta_lo)

      # if it is the last iteration, mark the intersection line
      if  iter == iter_limit
         plot!([x_lo,x_up],[y_lo,y_up],legend = false)
      end

      # update the new camber line interpolation data

        #examinexy[i,1:4] = [x_up,y_up,x_lo,y_lo]
        ad.c_x[i+1] = (x_up + x_lo)/2
        ad.c_y[i+1] = (y_up + y_lo)/2
      end
      # end of 1 whoel iteration to create a line camberline
      # print(slp_1[ad.n_intp + 1])

      # juedgement for the convergence :
      # average change of distance of the iteration on single camber point
      # to the chord length smaller
      chord  = sqrt((ad.c_x[1]-ad.c_x[end])^2 + (ad.c_y[1]-ad.c_y[end])^2)
      ds = 0
       for i = 1: ad.n_intp+2
          ds = ds + sqrt(((ad.c_x[i]-c_x_old[i])^2 + (ad.c_y[i]-c_y_old[i])^2)/ad.n_intp)
       end
      if ds/chord <= 1e-5
         print("camberline iter finished \n")
         break # break itertaion loop
      end

    end
   return ad
end
function num_camber3(cd::circle_data,ad::airfoil_data,iter_limit::Int64)

   # ===== Start of initializationin the method =========
   # initialize camber point data set in ad
   ad.c_x = zeros(Float64,ad.n_intp+2)
   ad.c_y = zeros(Float64,ad.n_intp+2)

    # starting initializing the slope of numerical camber line

    # initialize the theta matrix
    numtheta_up = zeros(Float64,ad.n_intp)
    numtheta_lo = zeros(Float64,ad.n_intp)
    numtheta_TE = 0.0
    numtheta_LE = 0.0

    # predefine function of coordinate needed in the method
    xi(theta)  = cd.xi_c +cd.k*cd.a*cos(theta)
    eta(theta) = cd.eta_c +cd.k*cd.a*sin(theta)

    r(theta)   = sqrt(xi(theta)^2 + eta(theta)^2)
    dr(theta)  = (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))

    dr2(theta) = (-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2/r(theta)^3

    dr3(theta) = (2*cd.xi_c*cd.k*cd.a*sin(theta) - 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^3-
                 1/2*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*cos(theta)-2*cd.eta_c*cd.k*cd.a*sin(theta))/r(theta)^3+
                 3/8*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^5


    x(theta)  = xi(theta)*(1 + 1/r(theta)^2)
    dx(theta) = -1*cd.k*cd.a*sin(theta)*(1 + 1/r(theta)^2) +
                xi(theta)*(-2/r(theta)^3)*dr(theta)
    dx2(theta) = xi(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*sin(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*cos(theta)*(1+1/r(theta)^2)
    dx3(theta) = xi(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*sin(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 + 1 + 1/r(theta)^2)+
                 cd.k*cd.a*cos(theta)*(6*dr(theta)/r(theta)^3)


    y(theta)  = eta(theta)*(1 - 1/r(theta)^2)
    dy(theta) = cd.k*cd.a*cos(theta)*(1 - 1/r(theta)^2) +
                eta(theta)*(2/r(theta)^3)*dr(theta)
    dy2(theta) = -1*eta(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*cos(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*sin(theta)*(1-1/r(theta)^2)
    dy3(theta) = -1*eta(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*cos(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 - 1 + 1/r(theta)^2)-
                 cd.k*cd.a*sin(theta)*(6*dr(theta)/r(theta)^3)

    # =================== curvature(kappa) module ==========================
    kappa(theta) = (dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)
    dkappa(theta) = (dx(theta)*dy3(theta)-dy(theta)*dx3(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)-
                    3/2*(dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(5/2)*(2*dx(theta)*dx2(theta)+2*dy(theta)*dy2(theta))

    total_point = 2*ad.n_intp+2
    kpa = zeros(Float64,total_point)
    phi = zeros(Float64,total_point)
    phi_edge = zeros(Float64,2)
    edge_index = 1
    for i = 1:total_point
        phi[i] = 2*pi/total_point*(i-1)
        kpa[i] = kappa(phi[i])
    end

    # find the local max point
    for i = 1:total_point
        L = i-1
        R = i+1
      if i == 1
        L = total_point
      elseif i == total_point
        R = 1
      end

      if kpa[i]-kpa[L]>=0 && kpa[i]-kpa[R]>=0
        phi_edge[edge_index] = phi[i]
        edge_index = edge_index + 1
        if edge_index == 3
            break
         end
       end
    end

    plot(phi,kpa,xlabel ="theta",ylabel = "/kappa",legend = false)
    savefig("kappa.png")

    print("kappa LE theta =",phi_edge[1],"\n")
    print("kappa TE theta =",phi_edge[2],"\n")
    # =========== end of kappa module ==================
    # redefine theta at LE/TE and update the point
    numtheta_TE = phi_edge[2]
    numtheta_LE = phi_edge[1]
    ad.c_x[1] = x(numtheta_TE)
    ad.c_y[1] = y(numtheta_TE)
    ad.c_x[ad.n_intp+2] = x(numtheta_LE)
    ad.c_y[ad.n_intp+2] = y(numtheta_LE)
    ad.x[1] = ad.c_x[1]
    ad.y[1] = ad.c_y[1]
    ad.x[ad.n_intp+2] = ad.c_x[ad.n_intp+2]
    ad.y[ad.n_intp+2] = ad.c_y[ad.n_intp+2]


    # use updated LE/TE position to redefine airfoil surface theta
    # these value will be used as initial guess for NEWTON Method
    lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
    up_range = 2*pi-lo_range
    for i = 1:ad.n_intp
     numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
     numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
     ad.x[i+1] = x(numtheta_up[i])
     ad.y[i+1] = y(numtheta_up[i])
     ad.x[end+1-i] = x(numtheta_lo[i])
     ad.y[end+1-i] = y(numtheta_lo[i])

    end
    plot(ad.x,ad.y,seriestype = :scatter,aspect_ratio =:equal,size = (1200, 800),dpi = 300)

    # ===== end of initializationin the method =========

    # predefine function needed for circle tangentil method :
    c(theta) = x(theta)
    d(theta) = -1*dx(theta)/dy(theta)
    dc(theta) = dx(theta)
    dd(theta) = -1*dx2(theta)/dy(theta)+dx(theta)*dy2(theta)/dy(theta)^2

    for i = 1:ad.n_intp

        print("camber point = ",i,"\n")
        print("upper point theta = ",numtheta_up[i],"\n")
        print("lower point theta = ",numtheta_lo[i],"\n")
        print("upper point x = ",x(numtheta_up[i]),"\n")
        print("lower point x = ",x(numtheta_lo[i]),"\n")
         a = -1*dy(numtheta_up[i])/sqrt(dx(numtheta_up[i])^2 + dy(numtheta_up[i])^2)
         b = x(numtheta_up[i])

        f(theta) =(2*a*(b-c(theta)))^2-4*(b-c(theta))^2*(a^2-1/(1+d(theta)^2))
        df(theta) = 2*(2*a*(b-c(theta)))*(-2*a*dc(theta))-
                    8*(b-c(theta))*(-1*dc(theta))*(a^2-1/(1+d(theta)^2))-
                    4*(b-c(theta))^2*(2*d(theta)*dd(theta))/(1+d(theta)^2)^2

        numtheta_lo[i] = mod(Newton_method(f,df,numtheta_lo[i]),2*pi)
        T = -1*(2*a*(b-c(numtheta_lo[i])))/(2*(a^2-1/(1+d(numtheta_lo[i])^2)))

        xc = b + a*T
        yc = y(numtheta_up[i])+dx(numtheta_up[i])/sqrt(dx(numtheta_up[i])^2 + dy(numtheta_up[i])^2)*T
        slp_up = dy(numtheta_up[i])/dx(numtheta_up[i])
        nom_up = (yc-y(numtheta_up[i]))/(xc-x(numtheta_up[i]))
        slp_lo = dy(numtheta_lo[i])/dx(numtheta_lo[i])
        nom_lo = (yc-y(numtheta_lo[i]))/(xc-x(numtheta_lo[i]))
        print("T = ",T,"\n")
        print("slp_up = ",slp_up,"\n")
        print("nom_up = ",nom_up,"\n")
        print("slp_lo = ",slp_lo,"\n")
        print("nom_lo = ",nom_lo,"\n")

         x_up = x(numtheta_up[i])
         y_up = y(numtheta_up[i])
         x_lo = x(numtheta_lo[i])
         y_lo = y(numtheta_lo[i])
         print("upper point x = ",x_up,"\n")
         print("lower point x = ",x_lo,"\n")
         ad.c_x[i+1] = (x_up + x_lo)/2
         ad.c_y[i+1] = (y_up + y_lo)/2
         plot!([x_lo,x_up],[y_lo,y_up],legend = false)
    end

   return ad
end
=#

function num_camber4(cd::circle_data,ad::airfoil_data)

    # ===== Start of initializationin the method =========
    # initialize camber point data set in ad
    ad.c_x = zeros(Float64,ad.n_intp+2)
    ad.c_y = zeros(Float64,ad.n_intp+2)
    ad.t   = zeros(Float64,ad.n_intp+2)

     # starting initializing the slope of numerical camber line

     # initialize the theta matrix
     numtheta_up = zeros(Float64,ad.n_intp)
     numtheta_lo = zeros(Float64,ad.n_intp)

     n_intpcan = 30*ad.n_intp
     numtheta_locan = zeros(Float64,n_intpcan)
     numtheta_TE = 0.0
     numtheta_LE = 0.0

     # predefine function of coordinate needed in the method
     xi(theta)  = cd.xi_c +cd.k*cd.a*cos(theta)
     eta(theta) = cd.eta_c +cd.k*cd.a*sin(theta)

     r(theta)   = sqrt(xi(theta)^2 + eta(theta)^2)
     dr(theta)  = (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))

     dr2(theta) = (-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))/(2*r(theta))-
                  1/4*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2/r(theta)^3

     dr3(theta) = (2*cd.xi_c*cd.k*cd.a*sin(theta) - 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))-
                  1/4*(-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))*
                  (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^3-
                  1/2*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))*
                  (-2*cd.xi_c*cd.k*cd.a*cos(theta)-2*cd.eta_c*cd.k*cd.a*sin(theta))/r(theta)^3+
                  3/8*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2*
                  (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^5


     x(theta)  = xi(theta)*(1 + 1/r(theta)^2)
     dx(theta) = -1*cd.k*cd.a*sin(theta)*(1 + 1/r(theta)^2) +
                 xi(theta)*(-2/r(theta)^3)*dr(theta)
     dx2(theta) = xi(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                  cd.k*cd.a*sin(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*cos(theta)*(1+1/r(theta)^2)
     dx3(theta) = xi(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                  cd.k*cd.a*sin(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 + 1 + 1/r(theta)^2)+
                  cd.k*cd.a*cos(theta)*(6*dr(theta)/r(theta)^3)


     y(theta)  = eta(theta)*(1 - 1/r(theta)^2)
     dy(theta) = cd.k*cd.a*cos(theta)*(1 - 1/r(theta)^2) +
                 eta(theta)*(2/r(theta)^3)*dr(theta)
     dy2(theta) = -1*eta(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                  cd.k*cd.a*cos(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*sin(theta)*(1-1/r(theta)^2)
     dy3(theta) = -1*eta(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                  cd.k*cd.a*cos(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 - 1 + 1/r(theta)^2)-
                  cd.k*cd.a*sin(theta)*(6*dr(theta)/r(theta)^3)

     # =================== curvature(kappa) module ==========================
     kappa(theta) = (dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)
     dkappa(theta) = (dx(theta)*dy3(theta)-dy(theta)*dx3(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)-
                     3/2*(dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(5/2)*(2*dx(theta)*dx2(theta)+2*dy(theta)*dy2(theta))

     total_point = 2*ad.n_intp+2
     kpa = zeros(Float64,total_point)
     phi = zeros(Float64,total_point)
     phi_edge = zeros(Float64,2)
     edge_index = 1
     for i = 1:total_point
         phi[i] = 2*pi/total_point*(i-1)
         kpa[i] = kappa(phi[i])
     end

     # find the local max point
     for i = 1:total_point
         L = i-1
         R = i+1
       if i == 1
         L = total_point
       elseif i == total_point
         R = 1
       end

       if kpa[i]-kpa[L]>=0 && kpa[i]-kpa[R]>=0
         phi_edge[edge_index] = phi[i]
         edge_index = edge_index + 1
         if edge_index == 3
             break
          end
        end
     end

     plot(phi,kpa,xlabel ="theta",ylabel = "/kappa",legend = false)
     savefig("kappa.png")

     print("kappa LE theta =",phi_edge[1],"\n")
     print("kappa TE theta =",phi_edge[2],"\n")
     # =========== end of kappa module ==================
     # redefine theta at LE/TE and update the point
     numtheta_TE = phi_edge[2]
     numtheta_LE = phi_edge[1]
     ad.c_x[1] = x(numtheta_TE)
     ad.c_y[1] = y(numtheta_TE)
     ad.c_x[ad.n_intp+2] = x(numtheta_LE)
     ad.c_y[ad.n_intp+2] = y(numtheta_LE)
     ad.x[1] = ad.c_x[1]
     ad.y[1] = ad.c_y[1]
     ad.x[ad.n_intp+2] = ad.c_x[ad.n_intp+2]
     ad.y[ad.n_intp+2] = ad.c_y[ad.n_intp+2]


     # use updated LE/TE position to redefine airfoil surface theta
     # these value will be used as initial guess for NEWTON Method
     lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
     up_range = 2*pi-lo_range
     for i = 1:ad.n_intp
      numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
      numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
      ad.x[i+1] = x(numtheta_up[i])
      ad.y[i+1] = y(numtheta_up[i])
      ad.x[end+1-i] = x(numtheta_lo[i])
      ad.y[end+1-i] = y(numtheta_lo[i])

     end
     plot(ad.x,ad.y,seriestype = :scatter,aspect_ratio =:equal,size = (1200, 800),dpi = 300)

     # predefined lower cancidate theta
     for i = 1:n_intpcan
       numtheta_locan[i] = mod(numtheta_TE- i*lo_range/(n_intpcan + 1), 2*pi)
     end

    for i = 1:ad.n_intp

       xu = x(numtheta_up[i])
       dxu = dx(numtheta_up[i])
       yu = y(numtheta_up[i])
       dyu = dy(numtheta_up[i])

       xl = zeros(Float64,n_intpcan)
       dxl = zeros(Float64,n_intpcan)
       yl = zeros(Float64,n_intpcan)
       dyl = zeros(Float64,n_intpcan)
       xc = zeros(Float64,n_intpcan)
       yc = zeros(Float64,n_intpcan)
       t1 = zeros(Float64,n_intpcan)
       t2 = zeros(Float64,n_intpcan)

       for j = 1:n_intpcan

           xl[j] = x(numtheta_locan[j])
           dxl[j] = dx(numtheta_locan[j])
           yl[j] = y(numtheta_locan[j])
           dyl[j] = dy(numtheta_locan[j])

           xc[j] = (xl[j]*dxl[j]/dyl[j]+yl[j]-yu-xu*dxu/dyu)/(dxl[j]/dyl[j]-dxu/dyu)
           yc[j] = -1*dxu/dyu*(xc[j]-xu) + yu

           t1[j] = sqrt((xu-xc[j])^2+(yu-yc[j])^2)
           t2[j] = sqrt((xl[j]-xc[j])^2+(yl[j]-yc[j])^2)

       end
       # find out the smallest different
       Tdif = 1
       T = 0
       index = 1
       for j = 1:n_intpcan
          tavg = (t1[j] + t2[j])/2
          if (t1[j]-t2[j])^2 <= Tdif
              Tdif = (t1[j]-t2[j])^2
              T = tavg
              index = j
              #print("t1 = ",t1[j],"t2 =" ,t2[j],"\n")
              #print("at $i point index = ", index,"\n")
          end
       end

       numtheta_lo[i] = numtheta_locan[index]
       x_up = x(numtheta_up[i])
       y_up = y(numtheta_up[i])
       x_lo = x(numtheta_lo[i])
       y_lo = y(numtheta_lo[i])

       ad.c_x[i+1] = (x_up + x_lo)/2
       ad.c_y[i+1] = (y_up + y_lo)/2
       plot!([x_lo,x_up],[y_lo,y_up],legend = false)

    end

    return ad

end
function num_camber5(cd::circle_data,ad::airfoil_data)

   # ===== Start of initializationin the method =========
   # initialize camber point data set in ad
   ad.c_x = zeros(Float64,ad.n_intp+2)
   ad.c_y = zeros(Float64,ad.n_intp+2)

    # starting initializing the slope of numerical camber line

    slp_1 = zeros(Float64,ad.n_intp+1) # camber sectional slope
    slp_2 = zeros(Float64,ad.n_intp)   # camber sectional average slope
    nom = zeros(Float64,ad.n_intp)     # camber sectional average normal
    c_x_old = zeros(Float64,ad.n_intp + 2)
    c_y_old = zeros(Float64,ad.n_intp + 2)

    # initialize the theta matrix
    numtheta_up = zeros(Float64,ad.n_intp)
    numtheta_lo = zeros(Float64,ad.n_intp)
    numtheta_TE = 0.0
    numtheta_LE = 0.0

    # predefine function of coordinate needed in the method
    xi(theta)  = cd.xi_c +cd.k*cd.a*cos(theta)
    eta(theta) = cd.eta_c +cd.k*cd.a*sin(theta)

    r(theta)   = sqrt(xi(theta)^2 + eta(theta)^2)
    dr(theta)  = (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))

    dr2(theta) = (-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2/r(theta)^3

    dr3(theta) = (2*cd.xi_c*cd.k*cd.a*sin(theta) - 2*cd.eta_c*cd.k*cd.a*cos(theta))/(2*r(theta))-
                 1/4*(-2*cd.xi_c*cd.k*cd.a*cos(theta) - 2*cd.eta_c*cd.k*cd.a*sin(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^3-
                 1/2*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))*
                 (-2*cd.xi_c*cd.k*cd.a*cos(theta)-2*cd.eta_c*cd.k*cd.a*sin(theta))/r(theta)^3+
                 3/8*(-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))^2*
                 (-2*cd.xi_c*cd.k*cd.a*sin(theta) + 2*cd.eta_c*cd.k*cd.a*cos(theta))/r(theta)^5


    x(theta)  = xi(theta)*(1 + 1/r(theta)^2)
    dx(theta) = -1*cd.k*cd.a*sin(theta)*(1 + 1/r(theta)^2) +
                xi(theta)*(-2/r(theta)^3)*dr(theta)
    dx2(theta) = xi(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*sin(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*cos(theta)*(1+1/r(theta)^2)
    dx3(theta) = xi(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*sin(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 + 1 + 1/r(theta)^2)+
                 cd.k*cd.a*cos(theta)*(6*dr(theta)/r(theta)^3)


    y(theta)  = eta(theta)*(1 - 1/r(theta)^2)
    dy(theta) = cd.k*cd.a*cos(theta)*(1 - 1/r(theta)^2) +
                eta(theta)*(2/r(theta)^3)*dr(theta)
    dy2(theta) = -1*eta(theta)*(6*dr(theta)^2/r(theta)^4-2*dr2(theta)/r(theta)^3) +
                 cd.k*cd.a*cos(theta)*(4*dr(theta)/r(theta)^3)-cd.k*cd.a*sin(theta)*(1-1/r(theta)^2)
    dy3(theta) = -1*eta(theta)*(18*dr(theta)*dr2(theta)/r(theta)^4-24*dr(theta)^3/r(theta)^5-2*dr3(theta)/r(theta)^3)+
                 cd.k*cd.a*cos(theta)*(6*dr2(theta)/r(theta)^3-18*dr(theta)^2/r(theta)^4 - 1 + 1/r(theta)^2)-
                 cd.k*cd.a*sin(theta)*(6*dr(theta)/r(theta)^3)

    # =================== curvature(kappa) module ==========================
    kappa(theta) = (dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)
    dkappa(theta) = (dx(theta)*dy3(theta)-dy(theta)*dx3(theta))/(dx(theta)^2+dy(theta)^2)^(3/2)-
                    3/2*(dx(theta)*dy2(theta)-dy(theta)*dx2(theta))/(dx(theta)^2+dy(theta)^2)^(5/2)*(2*dx(theta)*dx2(theta)+2*dy(theta)*dy2(theta))

    total_point = 2*ad.n_intp+2
    kpa = zeros(Float64,total_point)
    phi = zeros(Float64,total_point)
    phi_edge = zeros(Float64,2)
    edge_index = 1
    for i = 1:total_point
        phi[i] = 2*pi/total_point*(i-1)
        kpa[i] = kappa(phi[i])
    end

    # find the local max point
    for i = 1:total_point
        L = i-1
        R = i+1
      if i == 1
        L = total_point
      elseif i == total_point
        R = 1
      end

      if kpa[i]-kpa[L]>=0 && kpa[i]-kpa[R]>=0
        phi_edge[edge_index] = phi[i]
        edge_index = edge_index + 1
        if edge_index == 3
            break
         end
       end
    end

    #plot(phi,kpa,xlabel ="theta",ylabel = "/kappa",legend = false)
    #savefig("kappa.png")
    print("kappa LE theta =",phi_edge[1],"\n")
    print("kappa TE theta =",phi_edge[2],"\n")
    # =========== end of kappa module ==================
    # redefine data at LE/TE
    numtheta_TE = phi_edge[2]
    numtheta_LE = phi_edge[1]
    ad.c_x[1] = x(numtheta_TE)
    ad.c_y[1] = y(numtheta_TE)
    ad.c_x[ad.n_intp+2] = x(numtheta_LE)
    ad.c_y[ad.n_intp+2] = y(numtheta_LE)
    ad.x[1] = ad.c_x[1]
    ad.y[1] = ad.c_y[1]
    ad.x[ad.n_intp+2] = ad.c_x[ad.n_intp+2]
    ad.y[ad.n_intp+2] = ad.c_y[ad.n_intp+2]
    slp_1[1] = dy(numtheta_TE)/dx(numtheta_TE)
    slp_1[ad.n_intp] = dy(numtheta_LE)/dx(numtheta_LE)


    # use updated LE/TE position to redefine airfoil surface theta
    # these value will be used as initial guess for NEWTON Method
    lo_range = mod(numtheta_TE-numtheta_LE,2*pi)
    up_range = 2*pi-lo_range
    for i = 1:ad.n_intp
     numtheta_up[i] = mod(i*up_range/(ad.n_intp+1) + numtheta_TE,2*pi)
     numtheta_lo[i] = mod(numtheta_TE- i*lo_range/(ad.n_intp+1), 2*pi)
     ad.c_x[i+1] = (x(numtheta_up[i]) + x(numtheta_lo[i]))/2
     ad.c_y[i+1] = (y(numtheta_up[i]) + y(numtheta_lo[i]))/2
     ad.x[i+1] = x(numtheta_up[i])
     ad.y[i+1] = y(numtheta_up[i])
     ad.x[end+1-i] = x(numtheta_lo[i])
     ad.y[end+1-i] = y(numtheta_lo[i])
    end
    # ===== end of initializationin the method =========


    # generate the initial camber using kappa
    for i = 1:ad.n_intp
        # find the pairing for upper point
        k_theta_up = numtheta_up[i]
        k = kappa(k_theta_up)
        f(theta) = kappa(theta)-k
        df(theta) = dkappa(theta)
        k_theta_lo = mod(Newton_method(f,df,numtheta_lo[i]),2*pi)

        x_up = x(k_theta_up)
        y_up = y(k_theta_up)
        x_lo = x(k_theta_lo)
        y_lo = y(k_theta_lo)
        ad.c_x[i+1] = (x_up + x_lo)/2
        ad.c_y[i+1] = (y_up + y_lo)/2

        plot!([x_lo,x_up],[y_lo,y_up],legend = false)

    end

    return ad
end
airfoil = airfoil_data(num_intp)
airfoil = joukowsky(data,airfoil)
#examinexy = zeros(Float64,airfoil.n_intp,4)
#examineslope = zeros(Float64,airfoil.n_intp)
plot(airfoil.x,airfoil.y,aspect_ratio =:equal,size = (1200, 800),dpi = 300)

plot(airfoil.x,airfoil.y,seriestype = :scatter,aspect_ratio =:equal,size = (1200, 800),dpi = 300)
#airfoil = num_camber3(data,airfoil,3)
#plot(airfoil.x,airfoil.y,seriestype = :scatter,aspect_ratio =:equal,size = (1200, 800),dpi = 300)
#airfoil = num_camber2(data,airfoil,5)
airfoil = num_camber4(data,airfoil)
#airfoil = num_camber5(data,airfoil)
scatter!([airfoil.c_x[1],airfoil.c_x[end]],[airfoil.c_y[1],airfoil.c_y[end]],markershapes =:star5)
#scatter!(airfoil.c_x[end],airfoil.c_y[end],ma = :star5)
plot!([airfoil.c_x[1],airfoil.c_x[airfoil.n_intp+2]],[airfoil.c_y[1],airfoil.c_y[airfoil.n_intp+2]])
plot!(airfoil.c_x,airfoil.c_y,markershape = :cross, aspect_ratio =:equal,size = (1200, 800),dpi = 300)
#savefig("airfoil625verticle.png")
#plot(airfoil.c_x,airfoil.c_y,aspect_ratio =:equal,size = (1200, 800),dpi = 300)
#
# transformation function

# potential flow solution
# input variable :
p0 = 10000.0
rho = 1.74
U = 2.0
alfa = 15/180*pi
Gamma = -10

detail = 50
circle_xi = zeros(Float64,detail,detail)
circle_eta = zeros(Float64,detail,detail)
circle_u = zeros(Float64,detail,detail)
circle_v = zeros(Float64,detail,detail)
airfoil_x = zeros(Float64,detail,detail)
airfoil_y = zeros(Float64,detail,detail)
airfoil_u = zeros(Float64,detail,detail)
airfoil_v = zeros(Float64,detail,detail)
airfoil_p = zeros(Float64,detail,detail)
airfoil_cp = zeros(Float64,detail,detail)

function circle_field(cd::circle_data,num_r::Int64,num_theta::Int64)

     # creation of coordinate system for circle plane
         r_min = cd.a*cd.k
         r_max = 25*r_min

         for i = 1:num_r
             for j = 1:num_theta
             r = (i-1)/(num_r-1)*(r_max-r_min)+r_min
             theta = (j-1)/(num_theta-1)*2*pi
             circle_xi[i,j] = cd.xi_c + r*cos(theta)
             circle_eta[i,j] = cd.eta_c + r*sin(theta)
             circle_u[i,j] = U*cos(alfa) - U*(r_min/r)^2*cos(2*theta-alfa) -
                             Gamma/(2*pi*r)*sin(theta)
             circle_v[i,j] = U*sin(alfa) - U*(r_min/r)^2*sin(2*theta-alfa) +
                             Gamma/(2*pi*r)*cos(theta)
             end
         end
  return  circle_xi,circle_eta,circle_u,circle_v

end


function airfoil_field(cd::circle_data,num_r::Int64,num_theta::Int64)

    for i = 1:num_r
        for j = 1:num_theta

        zeta =  circle_xi[i,j] + circle_eta[i,j]*im
        z = zeta + (cd.xi_s + cd.eta_s*im)^2/zeta
        airfoil_x[i,j] = real(z)
        airfoil_y[i,j] = imag(z)

        transform = 1/(1-(cd.xi_s + cd.eta_s*im)^2/zeta^2)
        velo = (circle_u[i,j] - circle_v[i,j]*im)*transform

        airfoil_u[i,j] = real(velo)
        airfoil_v[i,j] = -1*imag(velo)
        end
    end
  return airfoil_x,airfoil_y,airfoil_u,airfoil_v
end

function steady_pressure(u::Array{Float64},v::Array{Float64},rho::Float64,p0::Float64)

     r_limit = size(airfoil_u)[1]
     theta_limit = size(airfoil_u)[2]
     for i = 1:r_limit
         for j  = 1:theta_limit
         airfoil_p[i,j] = (p0 + rho*U^2 -rho*(u[i,j]^2 + v[i,j]^2)/2)
         end
     end
     return airfoil_p
end
(circle_xi,circle_eta,circle_u,circle_v) = circle_field(data,detail,detail)
(airfoil_x,airfoil_y,airfoil_u,airfoil_v) = airfoil_field(data,detail,detail)
airfoil_p = steady_pressure(airfoil_u,airfoil_v,rho,p0)
plot(airfoil_x[1,:],(airfoil_p[1,:].-p0)/(1/2*rho*U^2),yflip = true)
#quiver(vec(circle_xi),vec(circle_eta),quiver = (vec(circle_u)*0.05,vec(circle_v)*0.05),
#       aspect_ratio = :equal,arrowsize = 0.01,linewidth = 0.1,size = (1200, 800),dpi = 300)
quiver(circle_xi[1,:],circle_eta[1,:],quiver = (circle_u[1,:]*0.05,circle_v[1,:]*0.05),
        aspect_ratio = :equal,arrowsize = 1,linewidth = 1.0,size = (1200, 800),dpi = 300)
#streamslice(vec(circle_xi),vec(circle_eta),vec(circle_u),vec(circle_v))
savefig("circlevelo.png")
plot(airfoil_x[1,:],airfoil_y[1,:])
quiver!(vec(airfoil_x[1:4,:]),vec(airfoil_y[1:4,:]),quiver = (vec(airfoil_u[1:4,:])*0.5,vec(airfoil_v[1:4,:])*0.5),aspect_ratio = :equal)
=#
