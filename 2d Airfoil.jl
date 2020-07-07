# Joulosky Airfoil Code
# Developed By Yi Tsung Lee (Joe)
# Version  V 1.2 (0703 2020)
# ==========================================================
# This code is capable of using joukowsky transform to generate 2d airfoil
# the main purpose is to generate 2d both edge rounded airfoil
# It will also calculate the corresponding 2D Potential Flow solution
# numcamber function has two version which use analytical joukowsky formulation
# to calculate the numerical camber line
# Or it is capable of providing TE/LE and surface point for the code to
# automatically calculate the numerical camber line
# ==========================================================

# ============== Pkg Used ==========================
#Pkg.add("Dierckx")
#using Plots
#gr()
using PyPlot
using Dierckx
using LinearAlgebra
using DelimitedFiles
# ==================================================
# ======= airfoil global variable definition ===============
k = 1.15     # Rounded Factor, >1.0 will generate both edge rounded airfoil
xi_s = 1.0  # Singular Point real value on Circle
eta_s = 0.0 # Singular Point imag value on Circle
xi_c = -0.3  # center Point real value on Circle : effect TE/LE curvature edge
eta_c = 0.5  # center Point real value on Circle : effect how much Camber
num_intp = 200 # define how much interpolation point using on a single surface
times = 70     # number of times for generating the candidate point to
               # calculate the camberline, 30 seems to work fine at this time
# ======= airfoil global variable definition ===============
# ======= flow global variable definition ===============
p0 = 10000.0  # free stream pressure
rho = 1.74
U = 2.0
alfa = 15/180*pi
Gamma = -10    # vortex strength

num_r = 200
num_theta = 200
field_scale = 10.0
# ======= flow global variable definition ===============
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
    numtheta ::Array{Float64}
    xi ::Array{Float64}
    eta::Array{Float64}
    theta_TE::Float64
    theta_LE::Float64
    numtheta_TE::Float64
    numtheta_LE::Float64
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
    data.numtheta = zeros(n_intp*2+2)
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

      # camber characteristic
      max_t::Float64
      max_tp::Float64 # max thickness position on camberline

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

figure(1,figsize = [12,8],dpi = 600)
subplot(111,aspect = "equal",xlabel = "ξ",ylabel = "η",title = "circls shape")
circle = circle_data(k,xi_s,eta_s,xi_c,eta_c)
circle = circle_generation(circle,num_intp)
plot(circle.xi,circle.eta,marker = ".")
plot([circle.xi[1],circle.xi_c],[circle.eta[1],circle.eta_c],marker = ".")
plot([circle.xi[num_intp+2],circle.xi_c],[circle.eta[num_intp+2],circle.eta_c],marker = ".")
savefig("circle.png")
PyPlot.display_figs()

figure(2,figsize = [12,8],dpi = 600)
subplot(111,aspect = "equal",xlabel = "x",ylabel = "y",title ="airfoil shape")
airfoil = airfoil_data(num_intp)
airfoil = joukowsky(circle,airfoil)
plot(airfoil.x,airfoil.y)
PyPlot.display_figs()

function num_camber(cd::circle_data,ad::airfoil_data,times::Int64)

    # ===== Start of initializationin the method =========
    # initialize camber point data set in ad
    ad.c_x = zeros(Float64,ad.n_intp+2)
    ad.c_y = zeros(Float64,ad.n_intp+2)
    ad.t   = zeros(Float64,ad.n_intp+2)

     # initialize the theta matrix
     numtheta_up = zeros(Float64,ad.n_intp)
     numtheta_lo = zeros(Float64,ad.n_intp)

     n_intpcan = times*ad.n_intp
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

     total_point = 2500    # this value assume 1000 is detailed enough
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

     figure(3,figsize = [12,8],dpi = 600)
     subplot(111,xlabel = "θ",ylabel = "κ",title = "Kappa distribution")
     plot(phi,kpa)
     PyPlot.display_figs()
     savefig("kappa.png")

     print("kappa LE theta =",phi_edge[1],"\n")
     print("kappa TE theta =",phi_edge[2],"\n")

     # =========== end of kappa module ==================
     # redefine theta at LE/TE and update the point
     numtheta_TE = phi_edge[2]
     numtheta_LE = phi_edge[1]
     cd.numtheta_TE = phi_edge[2]
     cd.numtheta_LE = phi_edge[1]
     cd.numtheta[1] = cd.numtheta_TE
     cd.numtheta[ad.n_intp + 2] = cd.numtheta_LE
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

      cd.numtheta[i+1] = numtheta_up[i]
      cd.numtheta[end+1-i] = numtheta_lo[i]

     end

     # predefined lower cancidate theta
     for i = 1:n_intpcan
       numtheta_locan[i] = mod(numtheta_TE- i*lo_range/(n_intpcan + 1), 2*pi)
     end

     figure(4,figsize = [12,8],dpi = 600)
     subplot(111,aspect = "equal",xlabel = "x",ylabel = "y")

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
          end

       end

       numtheta_lo[i] = numtheta_locan[index]
       x_up = x(numtheta_up[i])
       y_up = y(numtheta_up[i])
       x_lo = x(numtheta_lo[i])
       y_lo = y(numtheta_lo[i])

       ad.t[i+1]  = sqrt((x_up-x_lo)^2+(y_up-y_lo)^2)/2
       ad.c_x[i+1] = (x_up + x_lo)/2
       ad.c_y[i+1] = (y_up + y_lo)/2
       ad.x[i+1] = x_up
       ad.y[i+1] = y_up
       ad.x[end+1-i] = x_lo
       ad.y[end+1-i] = y_lo

       plot([x_lo,x_up],[y_lo,y_up])
    end
    plot(ad.c_x,ad.c_y)
    plot(ad.x,ad.y)
    title("camber pairing line")
    savefig("camber pairing line.png")
    PyPlot.display_figs()

    # =============== camber line characteristic ===========
    max_tp_index = 0
    max_t = 0
    for i = 1:ad.n_intp
      if ad.t[i+1] >= max_t
         max_t =  ad.t[i+1]
         max_tp_index = i+1
      end
    end

    ad.max_t = max_t
    ad.max_tp = ad.c_x[max_tp_index]

    # ======================================================
    # ================ examine module ======================
     err_angle = zeros(Float64,ad.n_intp)
     for i = 1:ad.n_intp

      c_slp1 = (ad.c_y[i+1]-ad.c_y[i])/(ad.c_x[i+1]-ad.c_x[i])
      c_slp2 = (ad.c_y[i+2]-ad.c_y[i+1])/(ad.c_x[i+2]-ad.c_x[i+1])
      c_nom = (ad.y[i+1]-ad.y[end+1-i])/(ad.x[i+1]-ad.x[end+1-i])
      c_slp = 1/2*(c_slp1 + c_slp2)
      err_angle[i] = abs((mod((atan(c_slp)-atan(c_nom)),pi)/pi*180-90))

     end
      err = (sum(err_angle))/ad.n_intp
      print("average error of camber degree=",err,"degree","\n")

      figure(5,figsize = [12,8],dpi = 600)
      subplot(111,aspect = "auto",xlabel = "x",ylabel = "error angle")
      plot(ad.c_x[2:end-1],err_angle)
      title("camber and pairing angle error distribution")
      savefig("error distribution.png")
      PyPlot.display_figs()


    # ================
    return ad

end
function num_general(x::Array{Float64},y::Array{Float64})

end

airfoil = num_camber(circle,airfoil,times)
figure(6,figsize = [12,8],dpi = 600)
subplot(111,aspect = "equal",xlabel = "x",ylabel = "y")
plot(airfoil.x,airfoil.y)
scatter([airfoil.c_x[1],airfoil.c_x[end]],[airfoil.c_y[1],airfoil.c_y[end]])
plot(airfoil.c_x,airfoil.c_y,marker = "x",color = "k")
plot([airfoil.c_x[1],airfoil.c_x[airfoil.n_intp+2]],[airfoil.c_y[1],airfoil.c_y[airfoil.n_intp+2]],color = "r")
title("airfoil plot")
legend(["airfoil","LE/TE point","camber line","chord line"])
savefig("airfoil.png")
PyPlot.display_figs()


# transformation function

# potential flow solution
circle_xi = zeros(Float64,num_r,num_theta)
circle_eta = zeros(Float64,num_r,num_theta)
circle_u = zeros(Float64,num_r,num_theta)
circle_v = zeros(Float64,num_r,num_theta)
airfoil_x = zeros(Float64,num_r,num_theta)
airfoil_y = zeros(Float64,num_r,num_theta)
airfoil_u = zeros(Float64,num_r,num_theta)
airfoil_v = zeros(Float64,num_r,num_theta)
airfoil_p = zeros(Float64,num_r,num_theta)
airfoil_cp = zeros(Float64,num_r,num_theta)

function circle_field(cd::circle_data,num_r::Int64,num_theta::Int64,dis::Float64)

     # creation of coordinate system for circle plane
         r_min = cd.a*cd.k
         r_max = dis*r_min

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
         airfoil_p[i,j] = (p0 + rho*U^2/2 -rho*(u[i,j]^2 + v[i,j]^2)/2)
         end
     end
     return airfoil_p
end

(circle_xi,circle_eta,circle_u,circle_v) = circle_field(circle,num_r,num_theta,field_scale)
(airfoil_x,airfoil_y,airfoil_u,airfoil_v) = airfoil_field(circle,num_r,num_theta)
airfoil_p = steady_pressure(airfoil_u,airfoil_v,rho,p0)


figure(7,figsize = [12,8],dpi = 600)
subplot(111,xlabel = "x", ylabel = "Cₚ")
plot(airfoil_x[1,:],(airfoil_p[1,:].-p0)/(1/2*rho*U^2))
axis([minimum(airfoil_x[1,:]),maximum(airfoil_x[1,:]),maximum((airfoil_p[1,:].-p0)/(1/2*rho*U^2)),minimum((airfoil_p[1,:].-p0)/(1/2*rho*U^2))])
title("static pressure coefficient distribution")
savefig("Static Cp distribution.png")
PyPlot.display_figs()

figure(8,figsize = [12,8],dpi = 600)
subplot(121,aspect = "equal",xlabel = "ξ", ylabel = "η")
#streamplot(airfoil_x)
quiver(circle_xi,circle_eta,circle_u,circle_v,linewidth = 0.001,width = 0.0004,headwidth = 5)

subplot(122,aspect = "equal",xlabel = "x", ylabel = "y")
plot(airfoil_x[1,:],airfoil_y[1,:],linewidth = 0.001)
quiver(airfoil_x,airfoil_y,airfoil_u,airfoil_v,linewidth = 0.001,width = 0.0004,headwidth = 5)
savefig("Flow Distribution.png")
PyPlot.display_figs()
