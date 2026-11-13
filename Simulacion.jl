using Plots
using DelimitedFiles

#Propiedades de los materiales
const k_acero=41.0 #Ternium
const pc_acero= 8131.0 * 434.0 #Ternium
const k_h2=0.183 #Tablas
const pc_h2=0.0899*14300.0 #Tablas

#Datos para la conductividad térmica radial del rollo
const k_oxido=40.0 #Inventado
const R_emulsion=6.67e-5 #Por datos
const emisividad_acero=0.08554 #Paper
const sigma_sb=5.67e-8 #Constante de Stefan-Boltzmann
const presion_contacto=8.0#Paper
const dureza_material=1133.86#Paper
const rugosidad_superficial=3.2215e-6#Paper
const tan_theta=0.1 #Inventado
const espesor_oxido=10e-6 #Paper
const varphi_contacto=presion_contacto/(dureza_material+presion_contacto) #Paper
const espesor_gap=42.7e-6*exp(-5e-2*presion_contacto)#Paper
const psi_cal=2 #Inventado

#Datos para el gas
const presion_gas=1.2e5 #Inventada
const masa_molar_h2=2.016e-3 #Tablas
const constante_gases=8.314 #Tablas
const u_gas = 11.76 #Ternium
const diametro_hidraulico=0.622 #Ternium
const longitud_canal=5.25 #Ternium
const coef_realce_radiacion=3.0# Coeficiente para Tg >= 800 K Por datos
const eta=2.0
const r_g = espesor_gap / ((1 - varphi_contacto) * k_h2)

#Modelo de propiedades del gas
const mu_h2 = 8.96e-6 #Tablas
rho_h2_t(T)=presion_gas*masa_molar_h2/(constante_gases*T)#Gases ideales

function h_conv_paper(tg_k)
    rho=rho_h2_t(tg_k)
    h=0.023*psi_cal*u_gas^0.8/diametro_hidraulico^0.2*k_h2*(mu_h2/rho)^(-0.8)*(1+(diametro_hidraulico/longitud_canal)^0.7)
    if tg_k>=800.0
        h*=coef_realce_radiacion
    else
        h=h
    end
    return h
end

#Calcular parametros para el rollo
function parametros_lambda_r(tg_k)
    S=a_metal*1e-2
    R_S2=S/(2*k_acero)
    R_D=(rugosidad_superficial/(1.13*k_acero*tan_theta))*varphi_contacto^(-0.94)
    if tg_k>=800.0
        return S,R_S2,espesor_oxido/k_oxido,R_D,S+espesor_gap+espesor_oxido,1.0
    else
        return S,R_S2,R_emulsion,R_D,S+espesor_gap,eta
    end
end


#Calcular lambda dependiente de la Temp
function calcular_lambda_r(tm,r_s2,r_cubierta,r_d,s_eff, etaprov)
    r_r=(2-emisividad_acero)/(4*(1-varphi_contacto)*emisividad_acero*sigma_sb*tm^3)
    r_paralelo=1/(1/r_r+1/r_g+1/r_d)
    return etaprov*s_eff/(2*r_s2+r_cubierta+r_paralelo)
end

#Actualizar toda la matriz de kr
function actualizar_kr!(kr,temp, idx_rollo,tg_k)
    r_s2,r_cubierta,r_d,s_eff,etaprov=parametros_lambda_r(tg_k)[2:end]
    @inbounds for idx in idx_rollo
        tm=temp[idx]+273.15
        kr[idx]=calcular_lambda_r(tm,r_s2,r_cubierta,r_d,s_eff,etaprov)
    end
end


#Hacer toda la geometria y clasificar flujos

function rectangulo()
    alturatotal=L*num_rollos+s_rollos*(num_rollos+1)+a_suelo+a_campint+a_atm+dv
    medioancho=R_ext+s_pared+a_campint+a_atm+dv
    anchdisc=UInt16(2*round(medioancho/dv))
    alturadisc=UInt16(round(alturatotal/dv))
    # aire/H2=0, rollo=1, campint=2, horno=3, suelo=4, placa=5
    material=zeros(alturadisc,anchdisc)
    material[:,end].=3
    material[:,1].=3
    material[1,:].=3

    num_indices_suelo=UInt16(round(a_suelo/dv))
    for i in 1:num_indices_suelo
        material[end-i+1,:].=4
    end

    num_indices_atm=UInt16(round(a_atm/dv))
    num_indices_campint=UInt16(round(a_campint/dv))

    for i in 1:num_indices_campint
        material[(2+num_indices_atm:end-num_indices_suelo),(end-num_indices_atm-i)].=2
        material[2+num_indices_atm:end-num_indices_suelo,num_indices_atm+i+1].=2
        material[1+num_indices_atm+i,(1+num_indices_atm+1:(end-1-num_indices_atm))].=2
    end
    
    H=material*-1
    num_indices_spared=UInt16(round(s_pared/dv))
    num_indices_srollos=UInt16(round(s_rollos/dv))
    num_indices_plate=UInt16(round(a_plate/dv))

    anchrollodisc=UInt16(round(R_ext/dv))
    alturarollodisc=UInt16(round(L/dv))

    rollo=zeros(alturarollodisc+num_indices_srollos,anchrollodisc)
    num_indices_rint=UInt16(round(R_int/dv))

    rollo[1:alturarollodisc,(1+num_indices_rint):end].=1
    rollo[alturarollodisc+1+UInt16(round((num_indices_srollos-num_indices_plate)/2)):alturarollodisc+UInt16(round((num_indices_srollos-num_indices_plate)/2))+num_indices_plate,(1+num_indices_rint):end].=5

    rollo_invertida=reverse(rollo,dims=2)
    rollo_t=[rollo_invertida rollo]
    rolloy,rollox=size(rollo_t)
    
    x_pos=2+num_indices_atm+num_indices_campint+num_indices_spared
    y_pos=2+num_indices_atm+num_indices_campint+num_indices_srollos

    for i in 1:num_rollos
        r=num_rollos-i+1
        rolloH=zeros(alturarollodisc+num_indices_srollos,anchrollodisc)
        rolloH[1:alturarollodisc,(1+num_indices_rint):end].=r
        rolloH[:,1+num_indices_rint].+=0.1
        rolloH[:,end].+=0.1
        rolloH[1,(1+num_indices_rint):end].+=0.2
        if r==1
            rolloH[alturarollodisc,(1+num_indices_rint):end].+=0.2
        else
            rolloH[alturarollodisc,(1+num_indices_rint):end].-=0.8
        end

        rolloH[alturarollodisc+1+UInt16(round((num_indices_srollos-num_indices_plate)/2)):alturarollodisc+UInt16(round((num_indices_srollos-num_indices_plate)/2))+num_indices_plate,(1+num_indices_rint):end].=-5

        rolloH_invertida=reverse(rolloH,dims=2)
        rolloH_t=[rolloH_invertida rolloH]

        H[y_pos:y_pos+rolloy-1,x_pos:x_pos+rollox-1].=rolloH_t
        material[y_pos:y_pos+rolloy-1,x_pos:x_pos+rollox-1].=rollo_t
        y_pos+=rolloy
    end

    return material,H
end

#Hacer la recurrencia de flujos laterales y superiores para cada rollo
function flujos(num_rollos,zeta=1/0.909)
    if num_rollos==1
        return [1],[1]
    end
    fd=zeros(num_rollos-1)
    fd[end]=1/(2+num_rollos)
    for k in num_rollos-2:-1:1
        fd[k]=zeta*fd[k+1]
    end
    fr_lateral=ones(num_rollos)
    for k in 2:num_rollos
        fr_lateral[k]=fr_lateral[k-1]*(1-fd[k-1])
    end
    fr_up=zeros(num_rollos)
    for k in 1:num_rollos-1
        fr_up[k]=fr_lateral[k]*fd[k]
    end
    fr_up[end]=fr_lateral[end]
    return fr_lateral,fr_up
end

#Definir a q parte del rollo le toca cada flujo
function H_suprema(fr_lat,fr_up,H,num_rollos)
    nf,nc=size(H)
    fz=ones(nf,nc)
    fr=ones(nf,nc)
    fz[H.<0].=1
    fr[H.<0].=1
    for r in num_rollos:-1:1
        fluz=fr_up[r]
        fluzr=fr_lat[r]
        fz[isapprox.(H,r+0.2,atol=1e-6).|isapprox.(H,r+0.3,atol=1e-6)].=fluz
        fr[isapprox.(H,r+0.1,atol=1e-6).|isapprox.(H,r+0.3,atol=1e-6)].=fluzr
    end
    return fr,fz
end


#Condicion frontera

function condfrontera!(temp_actual,temphorno)
    temp_actual[1,:].=temphorno
    temp_actual[:,1].=temphorno
    temp_actual[:,end].=temphorno
    temp_actual[end,:].=temphorno
end


#Aqui viene diff finitas
function paso!(tn,tnp1,es_algocond,kr,kz,pc,r_rec,h_gas,tg,dt)
    nf,nc=size(tn)
    @inbounds for j in 2:nc-1
        r=abs(r_rec[j])
        crecimiento=sign(r_rec[j])
        jmayor=UInt16(j+crecimiento)
        jmenor=UInt16(j-crecimiento)
        fac_mayor=(r+dv_m/2)/r
        fac_menor=(r-dv_m/2)/r

        for i in 2:nf-1
            if es_algocond[i,j]==0
                continue
            end
            
            temp_act=tn[i,j]
            q=0.0

            if es_algocond[i,jmayor]==0
                q+=fr[i,j]*h_gas*fac_mayor*(tg-temp_act)/dv_m
            else
                q+=kr[i,j]*fac_mayor*(tn[i,jmayor]-temp_act)/dv_m^2
            end

            if es_algocond[i,jmenor]==0
                q+=fr[i,j]*h_gas*fac_menor*(tg-temp_act)/dv_m
            else
                q+=kr[i,j]*fac_menor*(tn[i,jmenor]-temp_act)/dv_m^2
            end

            if es_algocond[i+1,j]==0
                q+=fz[i,j]*h_gas*(tg-temp_act)/dv_m
            else
                q+=kz[i,j]*(tn[i+1,j]-temp_act)/dv_m^2
            end

            if es_algocond[i-1,j]==0
                q+=fz[i,j]*h_gas*(tg-temp_act)/dv_m
            else
                q+=kz[i,j]*(tn[i-1,j]-temp_act)/dv_m^2
            end

            tnp1[i,j]=temp_act+dt*q/pc[i,j]
        end
    end
    return tnp1
end

#Temp del gas
datos = readdlm("Temperaturagas_A424_Tipo3.txt", '\t', Float64, skipstart=1)
const tiempos_tg = datos[:, 1] .* 3600.0
const temps_tg = datos[:, 2] .* 1.0
const temp_min = minimum(temps_tg)
const temp_max= maximum(temps_tg)
const temp_crit = temp_max
const t_max=maximum(tiempos_tg)

function t_horno(t)
    if t<=tiempos_tg[1]
        return temps_tg[1]
    end
    if t>=tiempos_tg[end]
        return temps_tg[end]
    end
    k=searchsortedlast(tiempos_tg,t)
    m=(temps_tg[k+1]-temps_tg[k])/(tiempos_tg[k+1]-tiempos_tg[k])
    return temps_tg[k]+m*(t-tiempos_tg[k])
end

#Variables iniciales que se pueden cambiar

const a_metal=0.15 #Ancho de metal de la lamina Ternium
const a_aire=0.00286 #Ancho de aire entre las laminas
const a_campint=12 #Ancho de campinterna
const a_atm=14 #Ancho de atm
const dv=2
const a_suelo=dv
const a_plate=32

const num_vueltas=430
num_rollos=4


const R_int=25.4
const L=125

const R_ext=R_int+num_vueltas*(a_metal+a_aire)-a_aire
println(R_ext - R_int)
const s_pared=24
const s_rollos=60


tiempo_total=3600*37.083 #Segundos
const tempsuelo=50

if a_plate>=s_rollos
    println("La placa no puede medir mas que la separacion de los rollos")
    exit()
end

#Resistencias iniciales
f_metal = a_metal / (a_metal + a_aire)
f_h2    = a_aire  / (a_metal + a_aire)

#Paralelo
kz_rollo = f_metal * k_acero + f_h2 * k_h2
pc_rollo = f_metal * pc_acero + f_h2 * pc_h2

propiedades = Dict(
    0 => (k_h2,     k_h2,     pc_h2),
    1 => (0.0,     kz_rollo, pc_rollo), # kr se sobrescribirá
    2 => (16.0,     16.0,     8000.0 * 500.0),
    3 => (25.0,     25.0,     2600.0 * 880.0),
    4 => (15.0,     15.0,     2000.0 * 880.0),
    5 => (16.0,     16.0,     8000.0 * 500.0)
)
function propiedades_grid(material,propiedades)
    nf,nc=size(material)
    kr=Array{Float64}(undef,nf,nc)
    kz=Array{Float64}(undef,nf,nc)
    pc=Array{Float64}(undef,nf,nc)
    for j in 1:nc,i in 1:nf
        kr[i,j],kz[i,j],pc[i,j]=propiedades[Int(material[i,j])]
    end
    return kr,kz,pc
end

material,H=rectangulo()

es_rollo=material.==1
es_placa=material.==5
es_campint=material.==2
es_algocond=es_rollo.|es_placa.|es_campint

#Inicializar variables
idx_rollo=findall(es_rollo)
nf,nc=size(material)
temp=fill(Float64(tempsuelo),nf,nc)
temp2=copy(temp)

kr,kz,pc=propiedades_grid(material,propiedades)
t_inicial_k = tempsuelo + 273.15 

actualizar_kr!(kr, temp, idx_rollo, t_inicial_k)

fr_lat,fr_up=flujos(num_rollos)
const fr,fz=H_suprema(fr_lat,fr_up,H,num_rollos)

const dv_m=dv*1e-2
centro_col=(nc+1)/2
r_rec=((1:nc).-centro_col).*dv_m

#Calcular paso de tiempo para que no se rompa
t_scan=300.0:10.0:1100.0
par_alta=parametros_lambda_r(1000.0)[2:end]
par_baja=parametros_lambda_r(300.0)[2:end]
h_max=0.01
kr_max=0.01
for t in t_scan
    kr_baja=calcular_lambda_r(t,par_alta...)
    kr_alta=calcular_lambda_r(t,par_baja...)
    h_alta=h_conv_paper(t)
    if h_alta>h_max
        global h_max=h_alta
    end
    if kr_baja>kr_max
        global kr_max=kr_baja
    end
    if kr_alta>kr_max
        global kr_max=kr_alta
    end
end
dt_max=pc_rollo/(2*(kr_max+kz_rollo)/dv_m^2+4*h_max/dv_m)
dt=0.9*dt_max
dt=6
#Donde evaluar cada cosa
n_upd=max(1,round(Int,30/dt))
n_pasos=max(1,round(Int,tiempo_total/dt))

#Calcular para los puntos calientes y frios de cada rollo
num_indices_atm=UInt16(round(a_atm/dv))
num_indices_campint=UInt16(round(a_campint/dv))
num_indices_spared=UInt16(round(s_pared/dv))
num_indices_srollos=UInt16(round(s_rollos/dv))

y_pos_rollo_ini=2+num_indices_atm+num_indices_campint+num_indices_srollos
rolloy=round(Int,L/dv)
anchrollodisc=round(Int,R_ext/dv)
num_indices_rint=round(Int,R_int/dv)
thickness=anchrollodisc-num_indices_rint
x_pos_ini=2+num_indices_atm+num_indices_campint+num_indices_spared

#Inicializar matrices por graficar
tempos = Vector{Float64}(undef, n_pasos)
temps_gas = Vector{Float64}(undef, n_pasos)
puntos_calientes = Array{Float64}(undef, num_rollos, n_pasos)
puntos_frios = Array{Float64}(undef, num_rollos, n_pasos)

yapaso=0

for paso_actual in 1:n_pasos
    t_h=t_horno((paso_actual-1)*dt)
    t_hk=t_h+273.15
    h_ref=h_conv_paper(t_hk)
    if (paso_actual-1)%n_upd==0
        actualizar_kr!(kr,temp,idx_rollo,t_hk)
    end
    paso!(temp,temp2,es_algocond,kr,kz,pc,r_rec,h_ref,t_h,dt)    
    global temp,temp2=temp2,temp
    print("\r",round(paso_actual/n_pasos*100,digits=1),"%")
    flush(stdout)


    y_pos_rollo=y_pos_rollo_ini

    for r in 1:num_rollos
        rr = num_rollos - r + 1
        y_mid=y_pos_rollo+div(rolloy-1,2)
        x_hot=x_pos_ini+2*anchrollodisc-1
        x_cold=x_hot-round(Int,(3/5)*thickness)

        puntos_calientes[rr, paso_actual] = temp[y_mid, x_hot]
        puntos_frios[rr, paso_actual] = temp[y_mid, x_cold]

        y_pos_rollo+=rolloy+num_indices_srollos
    end
    if minimum(puntos_frios[:, paso_actual])>=0.97*temp_max && yapaso==0
        global yapaso=(paso_actual-1)*dt
    end
    tempos[paso_actual] = (paso_actual-1) * dt
    temps_gas[paso_actual] = t_h
end

#Graficar heatmap
temp[.!(es_rollo.|es_placa.|es_campint)].=NaN
p1=heatmap(temp,yflip=true,color=:thermal,title="Temperatura (°C)",xlabel="X",ylabel="Y",colorbar=true)

nf,nc=size(material)
x_seg=Float64[]
y_seg=Float64[]

#Graficar lineas de contorno
for i in 1:nf
    for j in 1:nc-1
        if material[i,j]!=material[i,j+1]
            push!(x_seg,j+0.5);push!(x_seg,j+0.5);push!(x_seg,NaN)
            push!(y_seg,i-0.5);push!(y_seg,i+0.5);push!(y_seg,NaN)
        end
    end
end

for i in 1:nf-1
    for j in 1:nc
        if material[i,j]!=material[i+1,j]
            push!(x_seg,j-0.5);push!(x_seg,j+0.5);push!(x_seg,NaN)
            push!(y_seg,i+0.5);push!(y_seg,i+0.5);push!(y_seg,NaN)
        end
    end
end

append!(x_seg,[0.5,nc+0.5,nc+0.5,0.5,0.5,NaN])
append!(y_seg,[0.5,0.5,nf+0.5,nf+0.5,0.5,NaN])
#Graficar todo
Plots.plot!(p1,x_seg,y_seg,color=:black,linewidth=2.0,legend=false)
display(p1)

#Graficar puntos calientes y frios de cada rollo
p2=plot(tempos./3600,temps_gas,label="Temperatura del gas",color=:black,linestyle=:dash,linewidth=2,xlabel="Tiempo (h)",ylabel="Temperatura (°C)",title="Evolución Térmica: Punto Caliente vs. Punto Frío",legend=:right,size=(1200,750))
for i in 1:num_rollos
    plot!(p2,tempos./3600,puntos_calientes[i,:],label="Punto Caliente #$i",color=:red,linewidth=2)
    plot!(p2,tempos./3600,puntos_frios[i,:],label="Punto Frío #$i",color=:blue,linewidth=2)
end
plot!(p2,legend=:outertopright)
display(p2)

print("\nTiempo de llegada a 97% de la temperatura máxima: $(yapaso/3600) horas\n")
