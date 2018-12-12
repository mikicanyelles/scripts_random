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
### WT 12-LOX

p000t1 = gantt.Task(name='Abstraccions d\'hidrogen amb altres frames',
                    start=dt.date(2018, 12, 31),
                    duration=60,
                    percent_done=10,
                    resources=[rMCN],
                    )
p000.add_task(p000t1)



p.make_svg_for_tasks(filename='test_phd.svg',
                     today=dt.date(2018, 12, 12),
                     start=dt.date(2018, 12, 10),
                     end=dt.date(2018, 12, 31))
