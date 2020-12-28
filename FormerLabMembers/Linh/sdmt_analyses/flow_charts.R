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
  i [label = 'First sitting scores from participants\nwith at least 1 sitting']
  j [label = 'All sitting scores from participants\nwith at least 20 sittings ']
  i1 [label = 'Cross-sectional Data']
  j1 [label = 'Longitudinal Data']
  i11[label = 'HV v. MS comparison\n(Fig. 4A)']
  i13[label = 'If written SDMT from\nthe same visit are available']
  i12[label = 'If all clinical, written SDMT\nand MRI data from the same visit\nare available']
  i13a[label = 'Comparison of app\nand written SDMT performance\n(Fig. 4B-D)']
  i12a[label = 'Correlation and elastic\nnet analyses\n(Fig. 5 - 7)']
  j11[label = 'Practice effects\nanalyses\n(Fig. 8)']
  j12[label = 'ICC analyses\n(Fig. 9)']
  j13[label = 'Clinically meaningful\nchange analyses\nfor app SDMT\n(Fig. 11)']

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
  i->i1
  j->j1
  i1->{i11,i12,i13}
  i12-> i12a
  i13-> i13a
  j1-> {j11,j12,j13}

}")

grViz("digraph {

# setting graph statement
  
graph [overlap = true, fontsize = 14]

 node [shape = box,
        fontname = Helvetica]

a [label = 'Clinic data (written SDMT, EDSS, NeurEx, CombiWISE)']
b [label = 'Calculate differences in days between each visit.\nCalculate differences in written SDMT\nand clinical outcome scores between each visit.']
c [label = 'Keep written SDMT and clinical outcome\ndifferences for visits that were\nless than or equal to 190 days.']
d [label = 'Select only only the first written SDMT, EDSS, NeurEx,\nand CombiWISE differences for analyses.\n(Fig. 10; Days differences mean: 165.9 days; Range: 7 - 190 days)']

a->b->c->d
}")
