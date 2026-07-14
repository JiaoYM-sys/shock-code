# ------------------------------  Model settings  ------------------------------
clear
shell               mkdir output ./output/chuli ./output/tension

# ------------------------------  Basic settings  ------------------------------
units                           metal
dimension                       3
boundary                        p p s
atom_style                      atomic
timestep                        0.001
neighbor                        1.0 bin
neigh_modify                    every 1 delay 0 check yes

variable                        Temp      equal 300
variable                        tstp      equal 0.001
variable                        tdamp     equal 100*${tstp}
variable                        pdamp     equal 1000*${tstp}
variable                        dump_every equal 1000
variable                        equi_step equal 50000

# --------------------------  Build geometry  -----------------------------
read_data                       model.data

region                          fei  block INF INF INF INF 806 INF
region                          ba   block INF INF INF INF INF 806
group                           fei  region fei
group                           ba   region ba

region                          free block INF INF INF INF INF 100
group                           free region free
group                           mobile subtract all fei

write_data                      model.lmp

# --------------------------  Potential  ----------------------------------
pair_style                      hybrid eam/alloy adp lj/cut 10
pair_coeff                      * * eam/alloy almg.liu.eam.alloy Al NULL Mg
pair_coeff                      * * adp Si_Au_Al.mod.txt        Al Si  NULL
pair_coeff                      2 3 lj/cut 0.06 2.68

# -------------------------  应力/速度云图计算  -----------------------------
compute                         p_energy mobile pe/atom
compute                         spa      mobile stress/atom NULL
compute                         volum    mobile voronoi/atom

variable                        spaxx atom c_spa[1]/(c_volum[1]+0.01)*1e-4
variable                        spayy atom c_spa[2]/(c_volum[1]+0.01)*1e-4
variable                        spazz atom c_spa[3]/(c_volum[1]+0.01)*1e-4
variable                        spaxy atom c_spa[4]/(c_volum[1]+0.01)*1e-4
variable                        spaxz atom c_spa[5]/(c_volum[1]+0.01)*1e-4
variable                        spayz atom c_spa[6]/(c_volum[1]+0.01)*1e-4

variable                        Mises  atom sqrt(((v_spaxx-v_spayy)^2+(v_spayy-v_spazz)^2+(v_spaxx-v_spazz)^2+6*(v_spaxy^2+v_spaxz^2+v_spayz^2))/2)
variable                        tresca  atom (v_spazz-0.5*(v_spaxx+v_spayy))*0.5

# --------- 1D 分箱统计 ---------
compute                         ch_id mobile chunk/atom bin/1d z lower 1 units box
fix                             ave_ch all ave/chunk 5 100 500 ch_id v_spazz v_Mises v_tresca file ./output/chuli/stressz.txt

compute                         ch_v mobile chunk/atom bin/1d z lower 1 units box
fix                             ave_v mobile ave/chunk 5 100 500 ch_v vz file ./output/chuli/velocity1.txt

# ---- 自由面速度 ----
variable                        free_z equal vcm(free,z)
fix                             free_out1 free print 100 "$(step) ${free_z}" file ./output/chuli/velocity2.txt screen no
# --------- 温度 ---------
compute         KE   all ke/atom 
compute                         Temp1 mobile temp
variable        KB   equal 8.625e-5 
variable        TEMP atom c_KE/1.5/${KB} 
# --------- 1D 分箱统计 ---------
compute ch_temp mobile chunk/atom bin/1d z lower 1 units box

fix ave_temp mobile ave/chunk 5 100 500 ch_id v_spazz v_TEMP file ./output/chuli/temp.txt
compute                         p    all   pressure Temp1
variable                        pz     equal c_p[3]
variable                        p_GPa  equal v_pz*1e-4        
fix                             data_out all print 100 "$(step) ${p_GPa} " file ./output/chuli/p.txt screen no title "step p_GPa"
# ----------------------------------------------------------------------- 跟踪---------------------------------------------------------------------------------
region                       genzong block INF INF INF INF 210 211
group                        genzong region genzong
fix                             ave_temp1 genzong ave/chunk 5 100 500 ch_id v_spazz v_TEMP file ./output/chuli/temp1.txt

# --------------------------------  冲击加载  -------------------------------
velocity                        fei set 0 0 0
fix                             1 all nve
fix                             2 fei  move linear 0 0 -6 units box

# --------------------------  数据采集变量  -------------------------------
region                          measure1 block INF INF INF INF 200 300
region                          measure2 block INF INF INF INF 500 600

compute                         ba_average mobile reduce ave vz             #原子平均速度（统计平均，含热运动）
variable                        up equal vcm(mobile,z)               # 质心速度（整体运动，热运动被抵消）

variable                        lz0_ba   equal 806         
variable                        V0_ba    equal lx*ly*v_lz0_ba
print                           "Initial ba-volume, V0_ba = ${V0_ba}"

compute                         vor_ba ba voronoi/atom
compute                         V_ba   ba reduce sum c_vor_ba[1]
variable                        V_ba_real equal c_V_ba        
variable                        v_ratio equal v_V_ba_real/v_V0_ba

# --------------------------  输出设置  -------------------------------
thermo                          100
thermo_style                    custom step c_Temp1 c_p v_V0_ba v_V_ba_real c_ba_average 

fix                             data_out ba print 100 "$(step) ${p_GPa} ${v_ratio} ${up} " file ./output/chuli/p-v-u.txt screen no title "step p_GPa v_ratio up"

dump                            result all  custom 1000 ./output/tension/dump.*.lammpstrj id type x y z vx vy vz v_TEMP v_spazz v_Mises c_p_energy v_tresca


# ---- 靶板 zz 应力（GPa） ----
compute                         stress_ba ba stress/atom NULL
compute                         pzz_ba    ba reduce sum c_stress_ba[3]
variable                        sigmazz equal c_pzz_ba/vol*1e-4

# ---- 层裂指标：最大拉应力 & 孔洞分数 ----
variable                        prin_min equal -(c_pzz_ba/v_V_ba_real*1e-4)

# ---- 主输出：应力/自由面速度/拉应力/孔洞分数 ----
fix                             out ba print 100 "$(step) ${sigmazz} ${prin_min}" file ./output/chuli/shock.txt screen no title "step sigmazz(GPa) prin_min(GPa)"

# --------------------------  运行  -------------------------------
run                             5000
unfix                           2
run                             50000