$cutoff =shift; #score 200 or 400 

# read in protein id <-> gene name 
$f=shift; #info file
open F,"<$f";
<F>;
while(<F>){
    ($id,$name)=(split/\s+/)[0,1];
    $id=substr($id,5);
    #print "$id\t$name\n";
    $i++;
    $h{$id}=$name
}
@k=keys %h;
#print "$i,",scalar(@k),"\n";
close F;

# read in interaction link file
<>;
while(<>){ #link file
   s/\s+$//;
   ($n1,$n2,$s)=(split/\s+/);
   if($s<$cutoff){next}
    $n1=substr($n1,5);
    $n2=substr($n2,5);
    if($n1 eq $n2){die("two names are the same:$n1,$n2\n")}
    print "$h{$n1}\t$h{$n2}\t$s\n";
}


