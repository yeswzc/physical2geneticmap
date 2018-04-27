#!/usr/bin/perl -w
use strict;

my ($geneticmap,$bim,$fasize);
die "perl $0 [geneticMap] [bim] [fa.sizes] >[bim]" unless @ARGV==3;
($geneticmap,$bim,$fasize)=@ARGV;
my %size;
open FS,"$fasize" or die $!;
while(<FS>){
        my @e=split;
        next if @e<2;
        $e[0]=~s/^[CHRchr0]+//g;
        $size{$e[0]}=$e[1];
}
close FS;
#should be sorted
my %gm;
open GM,"$geneticmap" or die $!;
<GM>;
while(<GM>){
        s/\r//g;
        chomp;
        my @e=split /,/,$_;
        $e[0]=~s/^[CHRchr0]+//g;
        push @{$gm{$e[0]}},[@e[1,3]];
#Chr1,2335202,1,0.277186787
#
}
close GM;

###
my ($currentChr1,$currentCM1,$currentCM1POS1);
my ($currentChr2,$currentCM2,$currentCM2POS2);
my $lastCM=0; #to check
open BIM,"$bim" or die $!;
open BIM,"$bim" or die $!;
while(<BIM>){
        my @e=split;
        if(! defined $currentCM2 or $e[0] ne $currentChr2){
                initia_cm($e[0]);
                $lastCM=0;
        }
        my $cm;
        if($e[3]<=$currentCM2POS2){
                $cm=basePoint($currentCM2POS2,$currentCM2,$e[3]);
        }elsif($e[3]>=$currentCM2POS2 && $e[3]<=$currentCM1POS1){
                $cm=inter($currentCM2POS2,$currentCM2,$currentCM1POS1,$currentCM1,$e[3]);
        }else{
                while($e[3]>=$currentCM1POS1 && @{$gm{$e[0]}} ){
                        ($currentChr2,$currentCM2,$currentCM2POS2)=($currentChr1,$currentCM1,$currentCM1POS1);
                        my $t=shift @{$gm{$e[0]}};
                        ($currentChr1,$currentCM1,$currentCM1POS1)=($e[0],$t->[1],$t->[0]);
                }
                if($e[3]>=$currentCM1POS1 && empty(@{$gm{$e[0]}}) && $e[3]>$currentCM1POS1){
                        ($currentChr2,$currentCM2,$currentCM2POS2)=($currentChr1,$currentCM1,$currentCM1POS1);
                        my $chrLen=$size{$e[0]};
                        my $chrCM=$chrLen*$currentCM1/$currentCM1POS1;
                        ($currentChr1,$currentCM1,$currentCM1POS1)=($e[0],$chrCM,$chrLen);
                        print STDERR "$e[0] reached the maximum\nThe chr$e[0] Len=$chrLen\nThe chr$e[0] CM=$chrCM\n";
                }
                if($e[3]>$currentCM1POS1){
                        print STDERR "Error! $e[3] > chromosome size($size{$e[0]})\n\n";exit;
                }
                $cm=basePoint($currentCM2POS2,$currentCM2,$e[3]);
        }
        if($cm<$lastCM){print STDERR "Error latteral position with larger CM",join("\t",@e),"\n"; last;}
        print "$e[0]\t$e[1]\t$cm\t$e[3]\t$e[4]\t$e[5]\n";
}
close BIM;
1;
exit;

sub empty{
        my @a=@_;
        if(@a){
                return 0;
        }else{
                return 1;
        }
}
sub initia_cm{
        my $chr=shift;
        my $t=shift @{$gm{$chr}};
        ($currentChr2,$currentCM2,$currentCM2POS2)=($chr,$t->[1],$t->[0]);
        if($currentCM2==0){
                $t=shift @{$gm{$chr}};
                ($currentChr2,$currentCM2,$currentCM2POS2)=($chr,$t->[1],$t->[0]);
        }
        $t=shift @{$gm{$chr}};
        ($currentChr1,$currentCM1,$currentCM1POS1)=($chr,$t->[1],$t->[0]);
}

sub basePoint{
        my ($x1,$y1,$x)=@_;
        my $y=$x*$y1/$x1;
        return($y);
        1;
}

sub inter{
        my ($x1,$y1,$x2,$y2,$x)=@_;
        my $x0=searchBase($x1,$y1,$x2,$y2);
        my $y=$y1*($x-$x0)/($x1-$x0);
        return $y;
        1;
}

sub searchBase{
        my ($x1,$y1,$x2,$y2)=@_;
        my $x0=($x2*$y1-$x1*$y2)/($y1-$y2);
        return($x0);
        1;
}
