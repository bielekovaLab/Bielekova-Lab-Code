library(DiagrammeR)

grViz("digraph {
  
  # setting graph statement
  
  graph [overlap = true, fontsize = 14]
  
  # setting node statements
  node [shape = box,
        fontname = Helvetica]
  a [label = 'All app SDMT data\n(154 MS, 39 HV)']
  b [label = 'Sittings with 2 trials']
  b1 [label = '90s trial']
  b2 [label = '75s trial']
  c [label = 'Sittings with 1 trial']
  c1 [label = '75s trial']
  c2 [label = '90s trial']
  d [label = 'Calculate minimum \nreliable testing time\n(Fig. 3)']
  e [label = 'Adjust the outlier\nscore between 2 trials\n(Fig. 2)']
  f [label = 'Average the scores\nbetween 2 trials']
  g [label = 'Make score comparable\nto a 90s trial&rsquo;s score:\n(total score/75) x 90']
  h [label = 'Obtain overall score for each sitting.']
  i [label = 'Cross-sectional data']
  j [label = 'Longitudinal data']

  a->{b,c}
  b->b1
  b->b2[color = 'red']
  b1->{d,e}
  b2->g[color = 'red']
  g->e[color = 'red']
  g->h[color = 'DodgerBlue']
  e->f->h
  c->c1 [color = 'DodgerBlue']
  c->c2 
  c1->g [color = 'DodgerBlue']
  c2->h 
  h->{i,j}

}")

grViz("digraph {

# setting graph statement
  
graph [overlap = true, fontsize = 14]

 node [shape = box,
        fontname = Helvetica]

a [label = 'Clinic data (written SDMT, EDSS, NeurEx, CombiWISE)']
b [label = 'Calculate differences in days between each visit.\nCalculate differences in outcome scores between each visit.']
c [label = 'Keep outcome differences for visits that were\nless than or equal to 190 days.']
d [label = 'Select only only the first written SDMT, EDSS, NeurEx,\nand CombiWISE differences for analyses.\n(Fig. 10; Days differences mean: 165.9 days; Range: 7 - 190 days)']

a->b->c->d
}")
