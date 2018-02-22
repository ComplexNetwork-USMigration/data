#!/usr/bin/perl -w


#Copyright (C) 2000,2001,2002 The Regents of the University of California.

#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Author: Priya Mahadevan (pmahadevan@cs.ucsd.edu)
#



use strict;
use warnings;

#stores the degree for each node
my @node_degree = ();

#subroutine that computes factorial.....
sub factorial {

  my $num = shift;

  my $prod = 1;
  while ($num > 1) {
    $prod *= $num;
    $num--;
  }
  return $prod;
}

#AS nos are not incremental.... so, assign incremental nos... but preserve the old
my %node_map = ();
my $newid = 0;

my @new_old_id = ();
my %degree_dist;

my %node_neighbors;

my $edge_count = 0;


#read from the given skitter file.... (just AS AS)

my $x_index;
my $y_index;

#we need our output in 4 different files...
# the first one has info related to #nodes, #edges, #avg node degree etc.
# the next file has info about avg nbr connectivity as a function of node
# degree.. This is used to plot the graph
# Next, we need to write out the cluster coeff of the graph and plot it...
# Finally, we write the node degree distribution to a file and plot this
# as well....

my @linput = @ARGV;

if ($#ARGV < 4) {
  print " usage cluster_metrics.pl input_file output_file nbr_out clus_out degree_out \n";
  exit;
}

open (IN_FILE, $ARGV[0]) or die "Could not open inputfile file \n";
open (OUT_FILE, ">$ARGV[1]") or die "Could not open output file file \n";
open (NBRCONN_FILE, ">$ARGV[2]") or die "Could not open output file for nbr conn\n";
open (CLUS_FILE, ">$ARGV[3]") or die "Could not open output file for clustering\n";
open (DEG_FILE, ">$ARGV[4]") or die "Could not open output file for degree\n";

open (RICH_FILE, ">rich_club") or die "Could not open rich club connectivity file \n";

open (DEGCCDF_FILE, ">deg.ccdf") or die "Could not deg ccdf  file \n";
open (NBRCCDF_FILE, ">nbr.ccdf") or die "Could not deg ccdf  file \n";
open (CLUSLOG_FILE, ">clus.log") or die "Could not open file for writing log of CC \n";

open (CLUS_CCDF, ">clus.ccdf") or die "Could not open file for writing ccdf of CC \n";

open (PK1K2_FILE, ">deg_deg_distr.log") or die "Could not open file for writing log of CC \n";

my $new_map_file = $ARGV[0] . ".new_map" ;

open (NEW_FILE, ">$new_map_file") or die "Could not open file for writing new mappings of AS numbers \n";

print NEW_FILE "#AS_NUMBER   NEW_ID \n";
while (<IN_FILE>)
{

  my @a = split();

  if ($a[0] eq $a[1]) {
    next;
  }

  $edge_count++;
  if (exists ($node_map{$a[0]})) {

    # it already exists !!
    $x_index = $node_map{$a[0]};
    $node_degree[$x_index]++;

  }
  else {
    #this node does not exist... add it
    %node_map = (%node_map, $a[0], $newid);
    print NEW_FILE "$a[0] $newid \n";

    #increment the degree of this node by 1.
    $node_degree[$newid]++;
    $new_old_id[$newid] = $a[0];
    $x_index = $newid;
    $newid++;
  }

  if (exists ($node_map{$a[1]})) {
    # it already exists !!
    $y_index = $node_map{$a[1]};
    $node_degree[$y_index]++;

  }
  else {

   #this node does not exist... add it
    %node_map = (%node_map, $a[1], $newid);
    print NEW_FILE "$a[1] $newid \n";

    #again increment node degree
    $node_degree[$newid]++;
    $new_old_id[$newid] = $a[1];
    $y_index = $newid;
    $newid++;
  }

  #Need to add the above values to node_neighbors hash !!
  $node_neighbors{$x_index}{$y_index} = $y_index;
  $node_neighbors{$y_index}{$x_index} = $x_index;

} #end of while....


# now print the degree of all the nodes.....
my $i;
my $node_count = $#node_degree + 1; #returns  size of the array

print OUT_FILE "Number of nodes = $node_count \n";
print OUT_FILE "Number of edges = $edge_count \n";


my $total = 0;
my $max_degree = 0;
my $deg_tilda = 0; 

my %degree_of_node ; # this is same as node_array, except that it is a hash !
# we need them both, thats why I've this redundancy !!

for ($i = 0; $i < $node_count; $i++) {

    $deg_tilda += ($node_degree[$i]) ** 2;
    $degree_of_node{$i} = $node_degree[$i];

    $total += $node_degree[$i];
    if ($node_degree[$i] > $max_degree) {
      $max_degree = $node_degree[$i];
    }

    $degree_dist{$node_degree[$i]} ++;

}

$deg_tilda = ($deg_tilda) / $total;
my $avg = $total / $node_count;

my $max_deg_ratio = $max_degree / ($node_count -1);
print OUT_FILE "Average node degree = $avg \n";
print OUT_FILE "Maximum node degree = $max_degree\n";
print OUT_FILE "Maximum degree ratio = $max_deg_ratio\n";

# compute rich-club connectivity... First, sort nodes based on dec degree
# the first "r" nodes form the members of the rich-club.
# rich-club connectivity = (# of links connecting rich club members ) / r * (r-1) / 2


my @rank = ();
my $r_index = 0;

foreach my $elem ( sort {$degree_of_node{$b} <=> $degree_of_node{$a} } keys %degree_of_node) {
  $rank[$r_index] = $elem; # stores the node id in dec order of degree....
  #print "$r_index $rank[$r_index] \n";
  $r_index++;
}

#For the first node in rank, rich club conn = 1;
print RICH_FILE "1 1\n" ;

my $rich_club_count = 0;
my $n_clique = 0; # max no of highest degree nodes forming a clique....

#rank goes from 1 to N, not 0 to n-1.....

# for every node in dec order of rank, compute rich club connectivity....
for ($r_index = 1; $r_index < $node_count; $r_index++) {

  my $n1 = $rank[$r_index];
  #my $rich_club_count = 0;

  # go down the rank list and see if edge exists....
  for (my $j = 0; $j < $r_index; $j++) {

    my $n2 =  $rank[$j];

    if (exists ($node_neighbors{$n1}{$n2}) ) {
      $rich_club_count++;
    }

  }
  # compute rich club count for every $n1... number of nodes in the club
  # is r_index + 1
  my $r_conn =  ($rich_club_count * 2) / ( ($r_index) * ($r_index + 1) );

  if ($r_conn == 1) {
    $n_clique = $r_index + 1; # node rank such that r-conn = 1
  }
  my $normalized_rank = ($r_index + 1 ) / $node_count;

  #print RICH_FILE "$r_index $normalized_rank $r_conn \n";
  printf RICH_FILE "%.5f %.6f \n", $normalized_rank, $r_conn ;
}

print OUT_FILE "Top clique size = $n_clique \n";


# We now have to compute degree-degree coorelation.... i.e number of nodes
# of degree k1 connected to nodes of degree k2.. This will be a matrix form, though
# we use hash tables.....

 
my %deg_deg_edges = (); # number of edges between nodes of degree k and k1
my $node_deg = 0;
my $nbr_deg = 0;

#also, find avg neighbor connectivity....
my %node_nbrs_total_degree;

foreach $i (keys %node_neighbors) {

 $node_deg = $node_degree[$i];

  foreach my $nbr (keys %{$node_neighbors{$i}} )   {

    $nbr_deg = $node_degree[$nbr];

    if (exists ($deg_deg_edges{$node_deg}{$nbr_deg}) ) {
      $deg_deg_edges{$node_deg}{$nbr_deg} =  $deg_deg_edges{$node_deg}{$nbr_deg} + 1;
   }
    else 
      {
	$deg_deg_edges{$node_deg}{$nbr_deg} = 1;
      }

    if (exists(	$node_nbrs_total_degree{$i}) ) {
      	$node_nbrs_total_degree{$i} += $nbr_deg;
    } else {
	$node_nbrs_total_degree{$i} = $nbr_deg;
    }

  }

}

my $tmp;
my %d_cor; # This is basically P(k)
my %d_d_cor; #P(k, k1) - both are defined in dorogovtsev's paper....

#now print deg_deg_edges of the nodes....

foreach $i ( sort {$a<=>$b} keys %deg_deg_edges) {

  foreach my $j (sort{$a<=$b} keys %{$deg_deg_edges{$i}} )   {
    $tmp += $deg_deg_edges{$i}{$j} / (2 * $edge_count);
    $d_d_cor{$i}{$j} = $deg_deg_edges{$i}{$j} / (2 * $edge_count);
    print PK1K2_FILE "$i $j $d_d_cor{$i}{$j} \n";
  }

  $d_cor{$i} = ($tmp * $avg) / $i;
  $tmp = 0;
  #print OUT_FILE "$i  $d_cor{$i} \n";

}


#this is temp ! I want a maximal random graph !!
#foreach $i ( sort {$a<=>$b} keys %deg_deg_edges) {
#  foreach my $j (sort{$a<=$b} keys %{$deg_deg_edges{$i}} )   {

#    my $mrg_prob = ($i*$d_cor{$i} * $j * $d_cor{$j}) / ($avg * $avg);
#    print "$i $j $d_cor{$i} $d_cor{$j} $mrg_prob \n";
#  }
#}


#print "Now printing node_deg_distribution P(k) \n";
# Also, we have to compute entropy of node degree distribution...
# defined as:  - Sum(x){ P(x)log2(P(x) }
# This will be stored as a value......

#Poisson aproximation to Binomial distribution....

# nu = (N - 1) * p , where N = number of nodes, p = ratio of number of links in the graph to the maximum possible links.. p = m / (n * (n - 1)) / 2


#Poisson approximation is given by P(k) = (nu) ^ k * e ^ (-nu)
#                                         ------------------
#                                                 k!

#Now, ideally k goes from 1 to N-1 (which is teh max degree).. however, the deno quickly tends to zero, so we can ignore P(k) for large values of k.

my $binomial_prob = (2 * $edge_count) / ($node_count * ($node_count - 1) );
my $nu = ($node_count - 1) * $binomial_prob;
my $bin_const = (2.71828) ** (-$nu);

my $binomial_entropy = 0;

my @binomial_distr;

# we only count to 50, the larger values are very small and can be ignored.
for (my $j = 1; $j < 50; $j++) {
  $binomial_distr[$j] = ( ($nu ** $j) * $bin_const) / (factorial($j) ) ;
#   print "$j    $distr[$j] \n";
 $binomial_entropy += ( ($binomial_distr[$j] ) * ( log($binomial_distr[$j]) / log(2) ) );
}

$binomial_entropy = -1 * $binomial_entropy;

print DEG_FILE "#deg  deg_distrn   log_deg log_deg_distrn \n";
my $entropy = 0;

my %p_tilda ; #defined as k * p(k) / avg degree
my $h_tilda =0;

 foreach my $j (sort {$a<=> $b} keys %d_cor )   {

   $entropy += ( ($d_cor{$j} ) * ( log($d_cor{$j}) / log(2) ) );
   #$entropy += ( ($d_cor{$j} ) * ( log $d_cor{$j} ) );

   $p_tilda{$j} = ($d_cor{$j} * $j ) / $avg;
 
   #print OUT_FILE "$j $p_tilda{$j} \n";
   if ($p_tilda{$j} != 0) {
       $h_tilda += ($p_tilda{$j} * (log($p_tilda{$j}) / log(2) ) );
     }

   print DEG_FILE "$j  $d_cor{$j}\t";

   printf DEG_FILE "%f %f \n", log($j), log($d_cor{$j});
 }

$h_tilda = -1 * $h_tilda;
 
# we now compute mutual Information !! 
#For that we first need h(k1,k2)

my $H = 0;

foreach $i (keys %d_d_cor) {

  foreach my $j (keys %{$d_d_cor{$i}} )   {
    $H += ($d_d_cor{$i}{$j} * (log($d_d_cor{$i}{$j}) / log(2) ) ) ;
  }
}

$H = -1 * $H;
my $max_MI = 2 * $h_tilda;

my $mutual_information = $max_MI - $H;

my $MI_ratio = $H / $max_MI;

#print ccdf file for node degrees.... 
my $ccdf = 0;
#open (TMP_FILE, ">temp_deg_chk.out") or die "Could not open inputfile file \n";
 foreach my $j (sort {$b<=> $a} keys %degree_dist )   {

   #print TMP_FILE "$j $degree_dist{$j} \n";
   $ccdf += ($degree_dist{$j} /  $node_count) ;
   print DEGCCDF_FILE "$ccdf  $j \n";
 }

$entropy = -1 * $entropy;

my $ent_ratio = $entropy / $binomial_entropy;
print OUT_FILE "entropy of node degree distribution = $entropy\n";
print OUT_FILE "Entropy of Binomial Node degree distribution (maximum) = $binomial_entropy\n";
print OUT_FILE "Ratio of entropy to maximum possible entropy = $ent_ratio\n";

print OUT_FILE "Mutual Information of degree degree distribution = $mutual_information \n" ;
print OUT_FILE "Max possible mutual_information = $max_MI \n\n";
print OUT_FILE "Mutual Information ratio = $MI_ratio \n";

# to find assortative coefficient !!
my $ass_sum = 0;
my $ass_sq = 0;
my $ass_prod = 0;

foreach $i (keys %node_neighbors) {
  
  my $src_deg = $node_degree[$i];

  foreach my $elem (keys %{$node_neighbors{$i}} )   {
    
    if ($elem < $i) {next;}
    my $dst_deg = $node_degree[$elem];

    $ass_sum +=  ( ($src_deg + $dst_deg) /2 );
    $ass_prod += ($src_deg * $dst_deg);
    $ass_sq += ( (($src_deg * $src_deg) + ($dst_deg * $dst_deg) )  / 2 );
  }

}

#print "ass_sum = $ass_sum \n";
#print "ass_prod = $ass_prod \n";
#print "ass_sq = $ass_sq \n";

my $ass_num = 0;
my $ass_den = 0;

$ass_num = ($ass_prod / $edge_count) - (($ass_sum / $edge_count)* ($ass_sum / $edge_count) ) ;
$ass_den = ($ass_sq / $edge_count) -  (($ass_sum / $edge_count)*($ass_sum / $edge_count) ) ;

my $ass_coeff = $ass_num / $ass_den;

#print "num = $ass_num \n";
#print "den = $ass_den \n";
print OUT_FILE "assortative coefficient for the graph = $ass_coeff\n";


# Clustering coefficient of the graph !!

my $avg_nbr_degree = 0;
my $deg;
my @tmp_nbrs;
my $nbr_cnt = 0;
my @clust_coeff;
my $cluster_coeff_of_graph = 0;

my $max_local_clustering = 0;

foreach $i (keys %node_neighbors) {

  
  $nbr_cnt = 0;

  foreach my $elem (keys %{$node_neighbors{$i}} )   {
    $tmp_nbrs[$nbr_cnt] = $elem;
    $nbr_cnt++;
  }

  #print "Node $i has $nbr_cnt neighbors ***** $node_degree[$i] \n";

  for (my $j = 0; $j < $nbr_cnt; $j++) {

    for (my $k = $j + 1; $k < $nbr_cnt; $k++) {

      my $src = $tmp_nbrs[$j];
      my $dst = $tmp_nbrs[$k];
      
      #print "$src $dst \n";

      # check if there is an edge between these two....
      if ( exists($node_neighbors{$src}{$dst})) {
	$deg++;
      }
    }
  }

  if ($nbr_cnt <= 1) {

    $clust_coeff[$i] = 0;
  }
  else {
    $clust_coeff[$i] = (2 * $deg) / ( ($nbr_cnt) * ($nbr_cnt - 1) ) ;
  }

  if ($max_local_clustering < $clust_coeff[$i]) {
    $max_local_clustering = $clust_coeff[$i];
  }

  print CLUS_CCDF "$clust_coeff[$i] \n";


    $deg = 0;
    $cluster_coeff_of_graph += $clust_coeff[$i];

} # for all nodes computed clustering coefficient.....




$cluster_coeff_of_graph =  $cluster_coeff_of_graph / $node_count ;

#print OUT_FILE "Mean clustering (previous) = $cluster_coeff_of_graph \n";

print OUT_FILE "Maximum local clustering = $max_local_clustering\n";

# We now computed weighted average of the cluster coefficient....

my $weighted_cluster_coeff = 0;
my $num = 0;
my $den = 0;

for ($i = 0; $i < $#clust_coeff; $i++) {

  $num += ($node_degree[$i] * ($node_degree[$i] - 1) * $clust_coeff[$i] );
  $den += ($node_degree[$i] * ($node_degree[$i] - 1)  );

}

#print "num = $num, den = $den i = $i \n";

$weighted_cluster_coeff = $num / $den;

print OUT_FILE "Weighted cluster coeff of the graph (weighted on degree) = $weighted_cluster_coeff \n";

#to plot the graph, we need to store C(k) as a function of degree and not
# nodes....

my %expt_clus_coeff = ();
$tmp = 0;

for ($i = 0; $i < $#clust_coeff; $i++) {

  $tmp = $node_degree[$i];

  if (!exists($expt_clus_coeff{$tmp}) ) {
    $expt_clus_coeff{$tmp} = $clust_coeff[$i];
  } 
  else {
    $expt_clus_coeff{$tmp} += $clust_coeff[$i];
  }

}

# ABove values must be averaged...divide each value by number of nodes with that degree....

my $mean_clustering = 0;
$cluster_coeff_of_graph = 0;
$max_local_clustering = 0;

foreach my $k (keys %expt_clus_coeff ) {

  if ($degree_dist{$k} != 0) {
    $expt_clus_coeff{$k} = $expt_clus_coeff{$k} / $degree_dist{$k};
  } else {
    $expt_clus_coeff{$k} = 0;
  }

  if ($max_local_clustering < $expt_clus_coeff{$k}) {
    $max_local_clustering = $expt_clus_coeff{$k};
  }

  $mean_clustering += ($expt_clus_coeff{$k} * $d_cor{$k});
  $cluster_coeff_of_graph += ($k * ($k-1) * $d_cor{$k} * $expt_clus_coeff{$k} );
}

#$cluster_coeff_of_graph = $cluster_coeff_of_graph / ;

#This is C_bar !!
print OUT_FILE "Mean clustering C_bar = $mean_clustering\n";

# we now find Cluster coeff using theoritical values.....
#  Fix k... (compute it for every k from 1 to max_degree..
#  for i, j from 2.... max_degree, 
#  Sum(i) Sum(j) (deg[i] - 1)*(deg[j] - 1)*P(i,j)*P(i,k)*P(k,j)
#                ----------------------------------------------
#                  deg[i]*deg[j]*P(i)*P(j)*P(k)
#


#my deg_deg_edges and degree_dist hashes in reality store Number and not Prob..
# to find prob, we have to divide by node_count

my $val = 0;
my %theo_clus_coeff = ();
my $theo_mean_clus = 0;
my $theo_mean_clus_coeff = 0;

$tmp = 0;

my $num_coeff;
my $den_coeff;
my $mean_sq_degree = 0;

my $t1;

foreach my $k (keys %degree_dist) {

   foreach $i (keys %degree_dist) {

    if ($i < 2) {
      next;
    }

    foreach my $j (keys %degree_dist) {

     if ($j < 2) {
       next;
     }

     if ( (! exists($deg_deg_edges{$i}{$j})) || (! exists($deg_deg_edges{$i}{$k}))
          || (! exists($deg_deg_edges{$j}{$k})) ) {
       next;
     }

     $num = ($i - 1)*($j - 1) * ($d_d_cor{$i}{$j}) * ($d_d_cor{$i}{$k}) * ($d_d_cor{$j}{$k});

     $den = ($i)*($j)*($d_cor{$j})*($d_cor{$i});

     if ($den != 0) {
       $val += ($num / $den);
       #$tmp += ($num_coeff / $den_coeff);
     }
     
  } #foreach $j

 } #foreach $i

 $mean_sq_degree += ( $k * $k * $d_cor{$k} );

  $theo_clus_coeff{$k} = ($val * $avg * $avg * $avg) / ($node_count * $k * $k * $d_cor{$k} * $d_cor{$k}) ;

   #$tmp += $val;

  $val = 0;
 
#print "for degree $k degree dist $degree_dist{$k} theoretical clus coeff = $theo_clus_coeff{$k} \n";

} #foreach $k

my $avg_theo_clus_coeff = 0;

foreach my $k (keys %theo_clus_coeff) {
  $theo_mean_clus += ($theo_clus_coeff{$k} * $d_cor{$k});
  $theo_mean_clus_coeff += ($k * ($k -1) * $theo_clus_coeff{$k} * $d_cor{$k}) / ($mean_sq_degree - $avg) ;

}

$cluster_coeff_of_graph = $cluster_coeff_of_graph / ($mean_sq_degree - $avg);




#mean_clustering is C_bar
#theo_mean_clustering is PKK for C_bar...
#we need the ratio of PKK for C_bar / C_bar !

#next, $cluster_coeff_of_graph is C
#Theor_mean_clus_coeff is PKK for C
# we need the ratio of PKK for C / C

#uncor_clus_coeff is PK for C_bar...
# we need PKK for C_bar  / PK for C_bar..... 



#print OUT_FILE "Theoretical Mean clustering = $theo_mean_clus\n";

my $ratio_theo_mean_clus_to_mean_clus = $theo_mean_clus / $mean_clustering;

print OUT_FILE "PKK / measured mean clustering ratio = $ratio_theo_mean_clus_to_mean_clus\n";

print OUT_FILE "Theoretical cluster coefficient = $theo_mean_clus_coeff\n";

#Also, find clust coeff for uncorrrelated networks...
#my $uncor_clus_coeff;

my $PK_C_bar = ( ($mean_sq_degree - $avg) ** 2) / ($node_count * $avg * $avg * $avg);

my $ratio_PKK_to_PK = $theo_mean_clus / $PK_C_bar;

print OUT_FILE "PKK / PK ratio = $ratio_PKK_to_PK \n";

my $theo_k_clus = $avg / ($node_count);

my $ratio_K_PK = $theo_k_clus / $PK_C_bar;

print OUT_FILE "K / PK ratio = $ratio_K_PK \n";

print OUT_FILE "Clustering coeff of the graph = $cluster_coeff_of_graph\n";

my $ratio_C_to_theo_C = $theo_mean_clus_coeff / $cluster_coeff_of_graph;

print OUT_FILE "measured to PKK ratio  = $ratio_C_to_theo_C \n";

print CLUS_FILE "#Deg ExptCC  PredictedCC UncorrelatedCC \n";


foreach my $k (sort {$a<=>$b} keys %degree_dist) {

  print CLUS_FILE "$k $expt_clus_coeff{$k} $theo_clus_coeff{$k}  $PK_C_bar $mean_clustering $theo_mean_clus $cluster_coeff_of_graph $theo_k_clus\n";

 if ( ($k != 0) && ($expt_clus_coeff{$k} != 0) && ($theo_clus_coeff{$k} != 0) ) {
   printf CLUSLOG_FILE "%f %f %f \n", log($k), log($expt_clus_coeff{$k}), log($theo_clus_coeff{$k});
 }
}


# Average neighbor connectivity....
# For each node of degree k1, find avg neighbor degree and then average this value
# for all nodes of degree k1.

# node_nbrs_total_degree already stores the sum of the degrees of all its neighbors...
# we need to average this (sum of nbr degree/# of neighbors) for each node
# and store it as a function of node degree....

my $nbr_degree = 0;
$avg_nbr_degree = 0;
$deg = 0;

my $ratio_for_avg_nbr_connectivity = 0;

# we don't want the total degree of all neighbors.. we need teh average.. so divide it by node degree.....
#also, divide avg neighbor connectivity by n - 1
foreach my $j (keys %node_nbrs_total_degree )   {
  $node_nbrs_total_degree{$j} = $node_nbrs_total_degree{$j} / ($node_degree[$j] * ($node_count -1) );
}


my %avg_nbr_connectivity;

foreach $deg (keys %degree_dist) {
  foreach $i (keys %node_nbrs_total_degree) {
    if ($deg == $node_degree[$i]) {
      if (exists($avg_nbr_connectivity{$deg}) ) {
	$avg_nbr_connectivity{$deg} += ($node_nbrs_total_degree{$i}); # / $deg);
      } else {
	$avg_nbr_connectivity{$deg} =  ($node_nbrs_total_degree{$i} ); # / $deg);
      }
    }
  } # for all nodes.....

} #for all degrees....

my $avg_avg_nbr_conn;
my $max_avg_nbr_conn = 0;
my $deg_for_max_conn = 0;

foreach $i (keys %avg_nbr_connectivity) {
  $avg_nbr_connectivity{$i} = $avg_nbr_connectivity{$i} / $degree_dist{$i};
  $avg_avg_nbr_conn += ($d_cor{$i} * $avg_nbr_connectivity{$i});

    if ($max_avg_nbr_conn < $avg_nbr_connectivity{$i} ) {
    $max_avg_nbr_conn =  $avg_nbr_connectivity{$i} ;
    $deg_for_max_conn = $i;
  }
}

my $avg_nbr_num = $mean_sq_degree / ($avg * ($node_count - 1));

$ratio_for_avg_nbr_connectivity = $avg_nbr_num / $avg_avg_nbr_conn; 

print OUT_FILE "Average average neighbor degree = $avg_avg_nbr_conn\n";
print OUT_FILE "Maximum average neighbor degree = $max_avg_nbr_conn\n";
print OUT_FILE "PK / measured avg neighbor degree ratio = $ratio_for_avg_nbr_connectivity\n";
print OUT_FILE "Degree for which avg nbr connectivity is maximum = $deg_for_max_conn\n";

#print "Now printing avg nbr conectivity \n";
#also, print out ccdf for avg nbr connectivity....

$ccdf = 0;
my $prev_val = -1;
my $tmp_cnt = 1;

# we don't want the total degree of all neighbors.. we need teh average.. so divide it by node degree.....
#foreach my $j (keys %node_nbrs_total_degree )   {
#  $node_nbrs_total_degree{$j} = $node_nbrs_total_degree{$j} / ($node_degree[$j);
#}


#this node_nbrs_total_degrees already stores the nbr connectivity for each node..


foreach my $j (sort {$node_nbrs_total_degree{$b}<=>$node_nbrs_total_degree{$a}}  keys %node_nbrs_total_degree )   {

  if ($prev_val == $node_nbrs_total_degree{$j}) {
    $tmp_cnt++;
  } else {
    $ccdf += ($tmp_cnt/$node_count) ;
    $tmp_cnt = 1;
    $prev_val = $node_nbrs_total_degree{$j};
    print NBRCCDF_FILE "$node_nbrs_total_degree{$j} $ccdf \n";
  }
   
 }


print NBRCONN_FILE "#degree avg_nbr_conn log_deg  log_avg_nbr_conn \n";
 foreach my $j (sort {$a<=> $b} keys %avg_nbr_connectivity )   {
    print NBRCONN_FILE "$j  $avg_nbr_connectivity{$j} \t";
    printf NBRCONN_FILE "%f %f \n", log($j), log($avg_nbr_connectivity{$j});
 }

