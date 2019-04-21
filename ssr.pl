#!/usr/bin/perl -w

# @author    Khaled Al-Shamaa <k.el-shamaa@cgiar.org>
# @version   2.0 <last update July 14, 2012>
# @copyright 2010-2019, ICARDA
# @license   GPLv3

# Get input FASTA file name as an argument from command line call
$filename = $ARGV[0];

# Script parameters
$min_base        = 2;
$max_base        = 6;
$min_repeat      = 2;
$max_split       = 3;
$max_comp        = 6;
$max_error_ratio = 0.10;

# Open FASTA file for read
unless(open(FASTAFILE,"<$filename")){
    print "Unable to open $filename: $!\n";
    exit 1; 
}

# Output header
push(@lines, "Filename\t$filename\n");
push(@lines, "Minimum base length\t$min_base\n");
push(@lines, "Maximum base length\t$max_base\n");
push(@lines, "Minimum repeats\t$min_repeat\n");
push(@lines, "Maximum split length\t$max_split\n");
push(@lines, "Maximum compound split\t$max_comp\n");
push(@lines, "Maximum error ratio\t$max_error_ratio\n");

push(@lines, "ID\tCategory\tBase\tStart\tEnd\tLength\tType\tError\tSSR\n");

# Read FASTA file contents
while(<FASTAFILE>){

    # Is it a header line?
    if(/>/){
    
        # If there is previous sequence then process it
        if ($seq) { 
           &find_ssr; 
        }
        
        # Read sequence name and clear the sequence string
        $id  = substr($_, 1, index($_, "\n")-1);
        chomp($_);
        $seq = '';

    } else {
    
        # Accumulate sequence string
        chomp($_);
        $seq .= $_;
        
    }
}

# Process the last sequence
&find_ssr;

close FASTAFILE;

# Print out results
foreach $line (@lines) {
	#print $line;
    if ($line =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t\%(.+)\t(.+)\n$/ ) {
        $id    = $1;
        $motif = $3;
        $start = $4;
        $end   = $5;

        $motif =~ s/[\[\]]//g;
        $motif =~ s/\|/ /g;

        print "$id\t$motif\t$start\t$end\n";
    }
}

# Clean end of script
exit 0;

sub find_ssr {
    $prev_position = 0;
    $prev_length   = 0;
    
    # The main regular expression which search for SSR in the given sequence
    while ($seq =~ m/(([ATCG]{$min_base,$max_base})(.{0,$max_split}\2){$min_repeat,})/g) {
    
        $position = length($`);
        $length   = length($&);

        $ssr  = $1;
        $base = $2;
        
        # Remove repeated nucleotides to find out error ratio
        $mutations   = $ssr;
        $mutations   =~ s/$base//g;
        $error_ratio = length($mutations)/$length;

        if ($error_ratio > 0) {
            $chk_seq = $ssr;
            
            for ($x=$min_base; $x<=$max_base; $x++) {
                if ($chk_seq =~ m/(([ATCG]{$x})(.{0,$max_split}\2){$min_repeat,})/) {
                    $chk_ssr    = $1;
                    $chk_base   = $2;
                    $chk_length = length($&);
                    
                    $mutations  = $ssr;
                    $mutations  =~ s/$chk_base//g;
                    $chk_ratio  = length($mutations)/$length;
                    
                    if ($chk_ratio < $error_ratio) {
                        $base        = $chk_base;
                        $ssr         = $chk_ssr;
                        $length      = $chk_length;
                        $error_ratio = $chk_ratio;
                    }
                }
            }
        }
                    
        # If the microsatellite motifs fit the regular expression of SSR but error 
		# ratio still exceed accepted threshold then try to cut one nucleotide from 
		# the candidate SSR and try to check it again
        while ($error_ratio > $max_error_ratio && length($ssr) > $min_repeat*$min_base) {
        
            $sub_seq = substr $ssr, 1;
            
            if ($sub_seq =~ m/(([ATCG]{$min_base,$max_base})(.{0,$max_split}\2){$min_repeat,})/){
            
                $position = $position + length($`);
                $length   = length($&);
                
                $ssr  = $1;
                $base = $2;

                $mutations = $ssr;
                $mutations =~ s/$base//g;
                
                $error_ratio = length($mutations)/$length;
            } else {
                $ssr = $sub_seq;
            }
        }
        
        if (length($base) == 6) {
            # Simplify hexa- repeats to di- or tri- if possible
            if ($base =~ m/(.{2})\1\1/) {
                $base = $1;
            } elsif ($base =~ m/(.{3})\1/) {
                $base = $1;
            }
        } elsif(length($base) == 4) {
            # Simplify tetra- repeats to di- if possible
            if ($base =~ m/(.{2})\1/) {
                $base = $1;
            }
        }

        # If the percentage of non-repeated nucleotides within accepted limits
        if ($error_ratio <= $max_error_ratio) {
            # Write down detected SSR information into output stack

            if ($ssr =~ m/^($base)+$/ ) { 
                $is_perfect = "Perfect"; 
            } else { 
                $is_perfect = "Imperfect"; 
            }
            
            $err = sprintf "%.2f", $error_ratio*100;
            $end = $position + $length;

            if (!($base =~ m/^([ATCG])\1{1,$max_base}$/)) {
                push(@lines, "$id\tSimple\t$base\t$position\t$end\t$length\t$is_perfect\t\%$err\t$ssr\n");
            
                # Check for possible compound SSR cases when new simple SSR start 
				# within accepted length of nucleotides between stretches of 
				# microsatellite
                if ($prev_position > 0 && $prev_position+$prev_length+$max_comp > $position) {
                    
                    $chk_max = length($base);
                    $is_same = 'no';
                    
                    for ($split=1; $split<=$chk_max; $split++){
                        $a = substr($base, 0, $split);
                        $b = substr($base, $split);
                        $c = $b . $a;
                        if ($c eq $prev_base) { 
                            $is_same = 'yes';
                        }
                    }

                    $total_length = $position + $length - $prev_position;
                    $comp_ssr     = substr $seq, $prev_position, $total_length;
                    $end          = $prev_position + $total_length;
                    
                    if ($is_same eq 'yes') {

                        $mutations   = $comp_ssr;
                        $mutations   =~ s/$prev_base//g;
                        $error_ratio = length($mutations)/$total_length;
                        $err         = sprintf "%.2f", $error_ratio*100;
                        
                        $is_perfect = "Imperfect";
                        $category   = "Simple";
                        $show_base  = $prev_base;
                        
                    } else {          
                        
                        $mutations   = $comp_ssr;
                        $mutations   =~ s/$base//g;
                        $mutations   =~ s/$prev_base//g;
                        $error_ratio = length($mutations)/$total_length;
                        $err         = sprintf "%.2f", $error_ratio*100;

                        if ($comp_ssr =~ m/^($prev_base)+($base)+$/ ) { 
                            $is_perfect = "Perfect"; 
                        } else { 
                            $is_perfect = "Imperfect"; 
                        }

                        $category  = "Compound";
                        $show_base = "[$prev_base|$base]";
                    }

                    # Remove the previous two simple SSR to be replaced by this new 
					# compound SSR
                    pop @lines;
                    
                    # Handle the issue of Compound in previous SSR
                    $temp = pop @lines;
                    if ($temp =~ m/\tCompound\t\[(.+)\]\t(.+)\t(.+)\t(.+)\t(.+)\t\%(.+)\t(.+)\n$/ ) {
                        # Append current base name to the compound list
                        $show_base         = '[' . $1 . '|' . $base . ']'; 
                        
                        # Start from beginning of previous compound
                        $prev_position     = $2;
                        
                        # Total length is current simple SSR end - previous compound 
						# SSR beginning
                        $total_length      = $end - $prev_position;
                        
                        # Actual number of mutations in previous compound SSR is: compound 
						# SSR length x error ration / 100
                        $prev_mutations    = $4 * $6 / 100;
                        
                        # Actual number of mutations in the appended simple SSR is: 
                        # length of appended stretch (end of current simple SSR - 
						# end of the previous compound SSR) * current simple SSR 
						# error ration / 100
                        $current_mutations = ($end - $3) * $err / 100;
                        
                        # Calculate new mutations ratio (i.e. total) by add the 
						# number of mutations in previous compound SSR to the 
						# mutations of current simple SSR including nucleotides 
						# between stretches then divided total by total length of 
						# new compound SSR
                        $mutations_ratio=($prev_mutations+$current_mutations)/$total_length;
                        
                        if ($mutations_ratio == 0) { 
                            $is_perfect = "Perfect"; 
                        } else { 
                            $is_perfect = "Imperfect"; 
                        }
                        
                        $err = sprintf "%.2f", $mutations_ratio*100;
                        $comp_ssr = substr $seq, $prev_position, $total_length; 
                    }
                    
                    push(@lines, "$id\t$category\t$show_base\t$prev_position\t$end\t$total_length\t$is_perfect\t\%$err\t$comp_ssr\n");
                }
                
                $prev_position = $position;
                $prev_length   = $length;
                $prev_base     = $base;
            }
        }
    }
}