==============================================================================
HOST-TRANSITION HISTORY  (branch-level, height-stratified)
==============================================================================
state annotation key : host_group
  builds/cdv-1200/20260919/clades/america2/seed12345/america2_host_group.trees: 10001 trees, 1000 discarded as burn-in
  builds/cdv-1200/20260919/clades/america2/seed54321/america2_host_group.trees: 10001 trees, 1000 discarded as burn-in
trees analysed       : 18002
transitions recorded : 566018
per tree             : mean 31.44, sd 6.40
tip-date check       : max |reconstructed - label| = 0.0000 yr (OK)

--- TIP CENSUS ---------------------------------------------------------------
  state             tips  % tips    oldest    newest
  procyonid           25   38.5%    1992.0    2019.3
  mustelid            14   21.5%    2004.5    2019.4
  felid               10   15.4%    1992.7    2016.0
  domestic_dog         9   13.8%    2004.5    2023.6
  wild_canid           7   10.8%    2013.2    2017.5

--- ROOT STATE (posterior) ---------------------------------------------------
  felid             52.1%  (9374)
  procyonid         16.8%  (3021)
  wild_canid        13.9%  (2503)
  domestic_dog      10.4%  (1877)
  mustelid           6.8%  (1227)

  felid holds the root in 52% of trees while carrying 15% of the tips.
  Oldest tips (within 1 yr of 1992.0): felid, procyonid
  READ: a minority state that also anchors the oldest tips is winning
  the root. Outbound transitions from it are inflated by that anchoring;
  treat its rank in the table as an upper bound, not an estimate.

--- TRANSITIONS, ranked  (cutoff for 'deep' = 2011.1) ------------
  transition                          n      %  median             IQR  %deep  %term
  felid -> procyonid              87103  15.4%  2006.5   1991.8-2013.3    64%    39%
  procyonid -> wild_canid         76113  13.4%  2013.0   2009.6-2015.8    28%    41%
  procyonid -> domestic_dog       75588  13.4%  2011.9   2003.6-2015.1    48%    54%
  procyonid -> mustelid           61058  10.8%  2013.8   2008.2-2016.3    35%    50%
  wild_canid -> domestic_dog      38925   6.9%  2012.8   2004.4-2017.1    39%    71%
  domestic_dog -> mustelid        38610   6.8%  2012.5   2003.7-2012.9    37%    28%
  wild_canid -> procyonid         34167   6.0%  2012.1   2005.5-2015.4    45%    47%
  mustelid -> procyonid           27601   4.9%  2010.9   2006.5-2013.5    52%    43%
  domestic_dog -> procyonid       24709   4.4%  2009.5   2003.8-2012.9    62%    32%
  felid -> domestic_dog           19145   3.4%  2005.9   1995.5-2009.8    76%    37%
  procyonid -> felid              16330   2.9%  2004.7   1991.7-2013.2    63%     1%
  felid -> mustelid               11444   2.0%  2002.9   1993.4-2006.2    90%     7%
  wild_canid -> mustelid          10942   1.9%  2009.1   1998.7-2012.7    57%    29%
  felid -> wild_canid              9990   1.8%  2003.1   1992.6-2008.2    82%     6%
  domestic_dog -> wild_canid       9872   1.7%  2010.0   2005.0-2014.3    53%    30%
  mustelid -> domestic_dog         7185   1.3%  2005.9   1998.3-2012.0    68%    52%
  wild_canid -> felid              5705   1.0%  1996.2   1988.9-2011.1    75%     3%
  domestic_dog -> felid            4838   0.9%  2001.3   1989.6-2009.0    83%     1%
  mustelid -> wild_canid           3538   0.6%  2005.9   1997.4-2009.5    83%     9%
  mustelid -> felid                3155   0.6%  1997.4   1989.4-2007.0    86%     1%
  %deep = share older than the cutoff. %term = share landing on a
  terminal branch, i.e. explaining one tip's label rather than any
  onward transmission.

--- FOCUS: felid -> procyonid ---------------------------------------------------
  n = 87103  (15.4% of all transitions)
  jump year   median 2006.5, IQR 1991.8-2013.3, range 1952.8-2018.3
  all others  median 2011.8, IQR 2003.6-2015.0
  share within 10 yr of its own oldest jump: 0%
  on terminal branches: 39%

  by decade (count, and share on terminal branches):
    1950s       1 (  0%) term   0%  #
    1960s      56 (  0%) term   0%  #
    1970s    1172 (  1%) term   0%  #
    1980s    7050 (  8%) term   1%  ###
    1990s   22879 ( 26%) term  53%  ###########
    2000s   21636 ( 25%) term   2%  ##########
    2010s   34309 ( 39%) term  61%  ################

  Two separated humps mean two processes are being pooled: deep
  jumps from root anchoring, and shallow ones from something else.
  If the shallow hump is mostly terminal, it is tip labels, not spread.

==============================================================================
Branch-level counts are a LOWER BOUND: a branch with two endpoints in
different states is scored once regardless of how many jumps the
continuous-time process actually made along it.
==============================================================================
