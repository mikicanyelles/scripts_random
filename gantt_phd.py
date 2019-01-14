# _*_ coding: utf-8 _*_


import datetime as dt
import gantt

# Font
gantt.define_font_attributes(fill='black',
                             stroke='black',
                             stroke_width=0,
                             font_family="Verdana")

# Festius
#gantt.add_vacations(dt.date(2018, 12, 25))
#gantt.add_vacations(dt.date(2018, 12, 26))
#gantt.add_vacations(dt.date(2019,  1,  1))

# Treballador
rMCN = gantt.Resource('Miquel')

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
p201 = gantt.Project(name='Distances module')
p202 = gantt.Project(name='Frame selector module')
p203 = gantt.Project(name='QM/MM adaptor module')


p00.add_task(p000); p00.add_task(p001)
p0.add_task(p00)
p20.add_task(p200); p20.add_task(p201); p20.add_task(p202); p20.add_task(p203)
p2.add_task(p20)
p.add_task(p0), p.add_task(p1), p.add_task(p2)


# MaR1
## 12-LOX + LTA4H
### WT 12-LOX (p000)


p000t0 = gantt.Task(name='H-abstraction by 12-LOX',
                    start=dt.date(2018, 2, 10),
                    duration=2000,
                    percent_done=10)

p000.add_task(p000t0)

p000t1 = gantt.Task(name='Epoxidation',
                    start=dt.date(2018, 11, 14),
                    duration=100,
                    percent_done=20)

p000.add_task(p000t1)

p000t2 = gantt.Task(name='Oxigenation',
                    start=dt.date(2019, 4, 1),
                    duration=180)

p000.add_task(p000t2)


### 2MT LTA4H (p001)

p001t0 = gantt.Task(name='Generation of mutants',
                    start=dt.date(2018, 12, 17),
                    duration=60)
                    #resources=[rMCN])
p001.add_task(p001t0)
    # Càlcul de la dinàmica i anàlisi d'aquesta. Es fa amb la intenció de generar
    # tants rotàmers com sigui possible de l'ARG382, per a després clusteritzar-los
    # i refer el docking amb el rotàmer més freqüent. A més a més, servirà per veure
    # com evoluciona el Zn i la seva esfera de coordinació.

p001t1 = gantt.Task(name='Study and validation of generated mutant',
                    start=dt.date(2019, 3, 1),
                    duration=365,
                    depends_of=[p001t0])
                    #depends_of=[p001t0])
                    #resources=[rMCN])
p001.add_task(p001t1)
    # Un cop conegut el rotàmer més freqüent durant la dinàmica, s'utilitza com a
    # estructura pel docking.










#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
p.make_svg_for_tasks(filename='test_phd.svg',
                     #today=dt.date(2018, 12, 14),
                     start=dt.date(2018, 2, 1),
                     end=dt.date(2020, 2, 1),
                     scale=gantt.DRAW_WITH_MONTHLY_SCALE)
