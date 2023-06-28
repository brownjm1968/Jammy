module Jammy

using DataFrames
using DataFramesMeta
using Statistics

#### constants ##################################
hbar = 6.582119569e-16  # [eV-s]
neutmass = 939.56542052e6 # [eV/c^2]
c = 299792458 # m/s
mperfermi = 1e-15
barnperfermi2 = 1e-30/1e-28

@doc """
## Description
Computes neutron cross section using the Reich-Moore formalism

## Arguments
* `egrid::Vector{Float64}`: energy grid in eV to compute cross sections 
* `respars::DataFrame`: resonance parameters DataFrame, with starting rows `Jpi`,`Elam`,`Gg`, and following rows other channels units of eV
* `targetspin`: quantum spin of target
* `targetmass`: target mass, relative to neutron mass
* `projspin`: quantum spin of projectile
* `projmass`: projectile mass, relative to neutron mass
* `ac`: the interaction radius in F

## Returns
A Float64 of level shift factor

## Examples
```julia
using Jammy 

egrid::Vector{Float64} = [0.01:0.0001:12.0;]
respars = DataFrame(
    Jpi = [4.0,   3.0],
    Elam= [4.2801, 10.34],
    Gg  = [0.0530,0.0550],
    Gn  = [0.0039,0.00466],
)

targetspin = 3.5
targetmass = 180.94803
projmass = 1 # neutron 
projspin = 0.5
ac = 8.1271  # fermi = 1e-15 m

# calc the cross section 
sigt,sign,sigx,sigc = Jammy.calcxs(egrid,respars,targetspin,targetmass,projspin,projmass,ac);
```
""" ->
function calcxs(egrid,respars,targetspin,targetmass,projspin,projmass,ac,bs_approx=false)
    # --- Default setup -----
    notBSA = 1 # B = S approximation
    if bs_approx 
        notBSA = 0
    end
    skip = 4
    numChan = ncol(respars) - skip
    numRes  = nrow(respars)

    # initialize arrays to zero
    sigx = zeros(size(egrid,1),1)
    sigt = zeros(size(egrid,1),1)
    sign = zeros(size(egrid,1),1)
    sigc = zeros(size(egrid,1),1)
    Rcc = zeros(numChan,numChan)im
    Xcc = zeros(numChan,numChan)
    Pc = zeros(numChan,numChan)
    Lc = zeros(numChan,numChan)im
    Sc = zeros(numChan,numChan)
    phic = zeros(numChan,numChan)
    sinphic2 = zeros(numChan,numChan) # sin(phi_c)^2
    sin2phic = zeros(numChan,numChan) # sin(2*phi_c)
    cos2phic = zeros(numChan,numChan) # cos(2*phi_c)
    # determine available J-pi
    availablejpi = []
    for i = 1:numRes 
        if i==1
            append!(availablejpi,respars.Jpi[i])
        elseif respars.Jpi[i] in availablejpi
            continue
        else
            append!(availablejpi,respars.Jpi[i])
        end
    end
    numJpi = size(availablejpi,1)

    println("Number of resonances:     ",numRes)
    println("Number of part. channels: ",numChan)
    println("Number of spin groups:    ",numJpi)
    ######### Energy loop ###########
    for iE = 1:size(egrid,1)
        # energy dependent
        k = wavenumber(projmass,targetmass,egrid[iE])
        rho = ktimesr(projmass,targetmass,egrid[iE],ac)
        # zero R-matrix

        ######### J-pi loop ###########
        for iJ = 1:numJpi
            jpi = availablejpi[iJ]
            resparsJpi = @chain respars begin
                @rsubset :Jpi == jpi
            end
            numResJpi = nrow(resparsJpi)
            J = abs(jpi)
            gJ = (2*J+1)/( (2*projspin+1)*(2*targetspin+1) )
            avGg = mean(resparsJpi.Gg)
            Rcc .= zero(Rcc[1,1])
            for iRes = 1:numResJpi
                L=resparsJpi.L[iRes]
                B=-L
                for iChan = 1:numChan
                    # Penetrability, shift factor, phase-shift
                    Pc[iChan,iChan] = penetrat(projmass,targetmass,egrid[iE],ac,L)
                    Sc[iChan,iChan] = shift(projmass,targetmass,egrid[iE],ac,L)
                    phic[iChan,iChan] = phaseshift(projmass,targetmass,egrid[iE],ac,L)
                    # log der

                    Lc[iChan,iChan] = (Sc[iChan,iChan] - B)*notBSA + 1im*Pc[iChan,iChan]
                    sinphic2[iChan,iChan] = sin(phic[iChan,iChan])^2
                    sin2phic[iChan,iChan] = sin(2*phic[iChan,iChan])
                    cos2phic[iChan,iChan] = cos(2*phic[iChan,iChan])
                    for jChan = 1:numChan
                        # R-matrix
                        iwidth,jwidth = 0,0
                        iwidth = sqrt(resparsJpi[iRes,iChan+skip]/2/penetrat(projmass,targetmass,resparsJpi.Elam[iRes],ac,L))
                        jwidth = sqrt(resparsJpi[iRes,jChan+skip]/2/penetrat(projmass,targetmass,resparsJpi.Elam[iRes],ac,L))
                        Rcc[iChan,jChan] += iwidth*jwidth/(resparsJpi.Elam[iRes] - egrid[iE] - (1im*avGg)/2.0)
                    end
                end
            end # level, channel-channel
            X = sqrt(Pc) * inv(Lc) * inv(inv(Lc)-Rcc) * Rcc * sqrt(Pc)
            # total
            for iChan = 1:numChan
                sigt[iE] += 4*pi/k^2 * gJ * (sinphic2[iChan,iChan] + imag.(X[iChan,iChan])*cos2phic[iChan,iChan] - real.(X[iChan,iChan])*sin2phic[iChan,iChan] )
                sign[iE] += 4*pi/k^2 * gJ * (sinphic2[iChan,iChan]*(1-2*imag.(X[iChan,iChan])) - real.(X[iChan,iChan])*sin2phic[iChan,iChan] + (sum(imag.(X)^2,dims=1)[iChan] + sum(real.(X)^2,dims=1)[iChan]))
                sigc[iE] += 4*pi/k^2 * gJ * (imag.(X[iChan,iChan]) - (sum(imag.(X)^2,dims=1)[iChan] + sum(real.(X)^2,dims=1)[iChan]))
            end
            # reaction
            sigx[iE] += pi/k^2 * gJ * sum(imag.(X)^2 + real.(X)^2)
        end # J-pi's 
        # convert to barns
        sigt[iE] *= barnperfermi2
        sign[iE] *= barnperfermi2
        sigx[iE] *= barnperfermi2
        sigc[iE] *= barnperfermi2
    end # energies
    return sigt,sign,sigx,sigc
end


@doc """
## Description
Computes the wavenumber (inverse of wavelength) in F

## Arguments
* `projmass`: projectile mass, relative to neutron mass
* `targetmass`: target mass, relative to neutron mass
* `energy`: the kinetic energy of projectile in eV

## Returns
A Float64 of wavenumber

## Examples
```julia
k = wavenumber(1.0,181.0,45.0)
```
""" ->
function wavenumber(projmass,targetmass,energy)
    k = 1/hbar*sqrt( 2*projmass*targetmass/(projmass+targetmass) * targetmass/(projmass+targetmass) * neutmass * energy )/c # [eV-s]^-1 * [eV-s/m] = [m]^-1
    k *= mperfermi # times [1/m] [m/F] = [1/F]
    return k
end

@doc """
## Description
Computes the unitless quantity of wavenumber times radius

## Arguments
* `projmass`: projectile mass, relative to neutron mass
* `targetmass`: target mass, relative to neutron mass
* `energy`: the kinetic energy of projectile in eV
* `ac`: the interaction radius in F

## Returns
A Float64 of wavenumber multiplied by radius

## Examples
```julia
rho = ktimesr(1.0,181.0,45.0,8.0)
```
""" ->
function ktimesr(projmass,targetmass,energy,ac)
    wavenumber(projmass,targetmass,energy)*ac
end

@doc """
## Description
Computes the penetrability for a given particle-pair 

## Arguments
* `projmass`: projectile mass, relative to neutron mass
* `targetmass`: target mass, relative to neutron mass
* `energy`: the kinetic energy of projectile in eV
* `ac`: the interaction radius in F
* `L`: the angular momentum (must be less than 5)

## Returns
A Float64 of penetrability

## Examples
```julia
P_c = penetrat(1.0,181.0,45.0,8.0,0)
```
""" ->
function penetrat(projmass,targetmass,energy,ac,L)
    p = ktimesr(projmass,targetmass,energy,ac)
    if L==0
        return p
    elseif L==1
        return p^3/(1 + p^2)
    elseif L==2
        return p^5/(9 + 3*p^2 + p^4)
    elseif L==3
        return p^7/(225 + 458*p^2 + 6*p^4 + p^6)
    elseif L==4
        return p^9/(11025 + 1575*p^2 + 135*p^4 + 10*p^6 + p^8)
    else
        error("L > 4 not allowed.")
    end
end

@doc """
## Description
Computes the phase shift factor for a given particle-pair 

## Arguments
* `projmass`: projectile mass, relative to neutron mass
* `targetmass`: target mass, relative to neutron mass
* `energy`: the kinetic energy of projectile in eV
* `ac`: the interaction radius in F
* `L`: the angular momentum (must be less than 5)

## Returns
A Float64 of phase shift factor

## Examples
```julia
phi_c = phaseshift(1.0,181.0,45.0,8.0,0)
```
""" ->
function phaseshift(projmass,targetmass,energy,ac,L)
    p = ktimesr(projmass,targetmass,energy,ac)
    if L==0
        return p
    elseif L==1
        return p - atan(p)
    elseif L==2
        return p - atan(3*p/(3-p^2))
    elseif L==3
        return p - atan(p*(15-p^2) / (15-6*p^2))
    elseif L==4
        return p - atan(p*(105 - 10*p^2)/(105 - 45*p^2 + p^4))
    else
        error("L > 4 not allowed.")
    end
end

@doc """
## Description
Computes the level shift factor for a given particle-pair 

## Arguments
* `projmass`: projectile mass, relative to neutron mass
* `targetmass`: target mass, relative to neutron mass
* `energy`: the kinetic energy of projectile in eV
* `ac`: the interaction radius in F
* `L`: the angular momentum (must be less than 5)

## Returns
A Float64 of level shift factor

## Examples
```julia
S_c = shift(1.0,181.0,45.0,8.0,0)
```
""" ->
function shift(projmass,targetmass,energy,ac,L)
    p = ktimesr(projmass,targetmass,energy,ac)
    if L==0
        return 0
    elseif L==1
        return -1/(1 + p^2)
    elseif L==2
        return -(18 + 3*p^2)/(9 + 3*p^2 + p^4)
    elseif L==3
        return -(675 + 90*p^2 + 6*p^4)/(225 + 458*p^2 + 6*p^4 + p^6)
    elseif L==4
        return -(44100 + 4725*p^2 + 270*p^4 + 10*p^6)/(11025 + 1575*p^2 + 135*p^4 + 10*p^6 + p^8)
    else
        error("L > 4 not allowed.")
    end
end

end # module Jammy
