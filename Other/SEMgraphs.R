library(DiagrammeR)
grViz("
digraph SEM {
graph [layout = neato,
overlap = true,
outputorder = edgesfirst]
node [shape = rectangle]

a [pos = '0,1!', label = 'x1']
b [pos = '2,1!', label = 'x2']
c [pos = '4,1!', label = 'x3']
d [pos = '6,1!', label = 'x4']
e [pos = '-1,3!', label = 'Xint', shape = ellipse]
f [pos='0,0!',label='Ex1', shape=ellipse] 
g [pos='2,0!',label='epsilon_x2', shape=ellipse] 
h [pos='4,0!',label='epsilon_x3', shape=ellipse] 
i [pos='6,0!',label='epsilon_x4', shape=ellipse] 
f->a[label='1']
g->b[label='1']
h->c[label='1']
i->d[label='1']
f->g [label = 'rho_xx']
e->a [label='1'] 
e->b [label='1'] 
e->c [label='1'] 
e->d [label='1'] 
}")

