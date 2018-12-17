# _*_ coding: utf-8 _*_


import datetime as dt
import gantt

# Font
gantt.define_font_attributes(fill='black',
                             stroke='black',
                             stroke_width=0,
                             font_family="Verdana")

# Festius
gantt.add_vacations(dt.date(2018, 12, 25))
gantt.add_vacations(dt.date(2018, 12, 26))
gantt.add_vacations(dt.date(2019,  1,  1))

# Treballador
rMCN = gantt.Resource('miki')

# Vacances/dies lliures personals
# rMCN.add_vacations(dfrom=datetime.date(2014, 12, 29),dto=datetime.date(2015, 1, 4))

# Projectes

p = gantt.Project(name='PhD')
p0 = gantt.Project(name='MaR1')
p00 = gantt.Project(name='12-LOX + LTA4H')
p000 = gantt.Project(name='WT 12-LOX')
p001 = gantt.Project(name='2MT LTA4H')
p1 = gantt.Project(name='12-LOX Patent')
p2 = gantt.Project(name='Tools')
p20 = gantt.Project(name='pyLOXr')
p200 = gantt.Project(name='RMSD module')

p00.add_task(p000); p00.add_task(p001)
p0.add_task(p00)
p20.add_task(p200)
p2.add_task(p20)
p.add_task(p0), p.add_task(p1), p.add_task(p2)


# MaR1
## 12-LOX + LTA4H
### WT 12-LOX (p000)

p000.add_task(gantt.Task(name='Abstraccions d\'hidrogen amb altres frames',
                    start=dt.date(2018, 12, 31),
                    duration=60,
                    percent_done=10,
                    resources=[rMCN]))


### 2MT LTA4H (p001)

p001t0 = gantt.Task(name='MD 2MT-LTA4H noEDH',
                    start=dt.date(2018, 12, 17),
                    duration=3,
                    percent_done=25)
                    #resources=[rMCN])
p001.add_task(p001t0)
    # Càlcul de la dinàmica i anàlisi d'aquesta. Es fa amb la intenció de generar 
    # tants rotàmers com sigui possible de l'ARG382, per a després clusteritzar-los 
    # i refer el docking amb el rotàmer més freqüent. A més a més, servirà per veure
    # com evoluciona el Zn i la seva esfera de coordinació.

p001t1 = gantt.Task(name='Redocking EDH a 2MT LTA4H',
                    start=dt.date(2018, 12, 20),
                    duration=4,                  
                    percent_done=0)
                    #depends_of=[p001t0])
                    #resources=[rMCN])
p001.add_task(p001t1)
    # Un cop conegut el rotàmer més freqüent durant la dinàmica, s'utilitza com a
    # estructura pel docking. 

p001t2 = gantt.Task(name='Docking mutant alternatiu',
                    start=dt.date(2018, 12, 17),
                    duration=2,
                    percent_done=0)
p001.add_task(p001t2)









#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
p.make_svg_for_tasks(filename='test_phd.svg',
                     today=dt.date(2018, 12, 14),
                     start=dt.date(2018, 12, 1),
                     end=dt.date(2018, 12, 31),
                     scale=gantt.DRAW_WITH_DAILY_SCALE)
