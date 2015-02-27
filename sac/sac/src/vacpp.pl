#!/usr/bin/perl -s
#############################################################################
#
# VAC PreProcessor (Copyleft) Gabor Toth 1994 
#
# Translates the dimension independent notation to Fortran 90 by expanding
# the Loop Annotation SYntax (LASY). Also translates COMMON, INCLUDE:, !SHIFT.
#
# Usage: make file.f                           # Called by the Makefile
#        vacpp.pl [-d=23] -                    # Interactive
#        vacpp.pl file.t > file.f              # Translation
#        vacpp.pl -maxlen=72 file.t > file.f   # Limit line length
#        vacpp.pl -O file.t > file.f           # Optimize !SHIFT regions
#############################################################################

$ndim=3; $ndir=3; 
$phi=-9; $z=-8;
$if_cd=1; $if_mc=0; $if_fct=0; $if_tvdlf=0; $if_tvd=0; 
$if_impl=0; $if_poisson=0; $if_ct=0; $if_gencoord=0; $if_resist=0; $if_rk=1;
$if_mpi=1;

# SETVAC READS UP TO THIS POINT

# Default maximum length for lines unless set by eg. "vacpp.pl -maxlen=72 ..."
$maxlen=78 unless $maxlen;

# For interactive use accept -d=12 flag
if($d){
    $d=~/([123])([123])/ || die "Incorrect -d flag value\n";
    $ndim=$1; $ndir=$2;
}

&definepatterns;

foreach $file (@ARGV) {
   &processfile($file, 'fh00');
}
#============================================================================
sub processfile {
   local($filename, $input) = @_;
   $input++;
   open($input, $filename) || die"Can't open $filename: $!\n";
   while (<$input>) {
      # PRINT EMPTY AND COMMENT LINES AS THEY ARE, EXCEPT FOR !SHIFT COMMENTS
      if(/^$/ || /^ *![^S]/){print; next}
      # HIDE QUOTED TEXT FROM TRANSLATION
      $_=&quotation($_);
      # PRINT TRAILING COMMENTS (commented out to protect !F77_ etc....)
      #if(s/ *(!.*)//){$comment=&unquote($1); /^ */; print "$&$comment\n"}
      # COLLECT CONTINUATION LINES
      next if &contline;
      # REPLACE TABS BY SPACES
      if(/\t/){
	  print STDERR "Warning: TABs are replaced by SPACEs\n" unless 
                                                               $tabwarn++;
          1 while s/\t+/' ' x (length($&)*8 - length($`)%8)/e;
      }
      # EXPAND THE LINE
      &processline($input);
      # PROCESS FILES INCLUDED BY "INCLUDE:filename"
      if (/INCLUDE:(.*)/) {&processfile($1,$input);next}
      # OPTIMIZE SHIFT STATEMENTS
      &shiftopt if $O && (/^ *!SHIFT/ || $shiftdef || $shiftproc);
      # PRINT THE LINE FORMATTED ACCORDING TO MAXLEN
      &printline if $_;
   }
   # Print common declarations collected from "COMMON,..." lines
   for $_ (values(%common)){&printline}; undef %common;
}
#===========================================================================
sub contline{
    # Join continuation lines into a single line. Uses global $contline
    if(/[^A-Z]\& *$/){$contline.=$_}
    else             {$_=$contline.$_;s/([^A-Z])\& *\n */$1/g;$contline=""}
}
#============================================================================
sub processline{
   local($input)=@_;
   # ATTACH A NUL-SEPARATOR CHARACTER TO THE BEGINNING OF LINE FOR LEFTMATCH
   $_=$spc.$_;
   # DO SUBSTITUTIONS FOR BLOCKS "{....}"
   while(/$ope/){
      &block($input);
      #TST print "Block:$_";
   }
   # DO SUBSTITUTION FOR EXPRESSIONS "...^pattern..."
   while(/$pat/){
      &expression;
      #TST print "Expr:$_";
   }
   # REMOVE NUL-SEPARATOR CHARACTERS
   s/$spc//g; 
   # REPLACE BREAK CHARACTERS BY NEWLINES AND INDENTATION
   if(/\\/){/^ */; $indent=$&; s/\\/\n$indent/g};
   # PROCESS "COMMON, ..." DECLARATIONS
   if(/^ *COMMON,/){&common};
}
#============================================================================
sub definepatterns{
   # Define special characters for the PreProcessor
   $ope='{'; $clo='}'; $bar='\|'; 
   $pat='\^'; $patid='[A-Z]'; $patchr='[A-Z&%]'; $uni='%';
   $brk='\\'; $sep='[,\+\-\*/:; ]'; $spc='~';
   $rbound=" ,;~)\n"; $lbound=" ,;~(\n";

   $nsubdefault=$ndim;
   # Define number of substitutes and subtitute strings for patterns
   # E.g. ^D -> 1,2 is defined by &patdef('D',2,'1','2','3')
   &patdef('LM'	,$ndim	,'ixGlo1:ixGhi1+2'	,'ixGlo2:ixGhi2+2','ixGlo3:ixGhi3+2'	);
   &patdef('SIDEADD'	,1	,'ixGlo^D:ixGhi^D+2'	);
   &patdef('SIDEADO'	,1	,'ixGlo^D:ixGhi^D+1'	);   
   &patdef('ND'	,1	,$ndim			);
   &patdef('NC'	,1	,$ndir			);
   &patdef('PHI',1	,$phi			);
   &patdef('Z'	,1	,$z			);
   if($phi>0){$pphi=$phi}else{$pphi=1};
   if($z>0){$zz=$z}else{$zz=1};
   &patdef('PPHI',1	,$pphi			);
   &patdef('ZZ'	,1	,$zz			);

   &patdef('DE'	,$ndim-1,2	,3		);
   &patdef('DE&',$ndim-1			);
   &patdef('DE%',$ndim-1,'^%2'	,'^%3'		);
   &patdef('SE'	,$ndim-1,'^LIM2:','^LIM3:'	);

   &patdef('D'	,$ndim	,1	,2	,3	);
   &patdef('D&'	,$ndim				);
   &patdef('DLOOP',$ndim       			);
   &patdef('DB'	,$ndim	,$ndim	,$ndim-1,$ndim-2);
   &patdef('DD'	,$ndim	,'^D' 	,'^D'	,'^D'	);
   &patdef('DD&',$ndim	,'^D&'	,'^D&'	,'^D&'	);
   &patdef('DDLOOP',$ndim,'^DLOOP','^DLOOP','^DLOOP');
   &patdef('D%'	,$ndim	,'^%1'	,'^%2'	,'^%3'	);
   &patdef('DL'	,$ndim	,'^LIM1','^LIM2','^LIM3');
   &patdef('S'	,$ndim	,'^LIM1:','^LIM2:','^LIM3:');
   &patdef('T'	,$ndim	,'^LLIM1:','^LLIM2:','^LLIM3:');
   &patdef('DLB',$ndim	,'^LIM'.$ndim,'^LIM'.($ndim-1),'^LIM'.($ndim-2));

   &patdef('CE'	,$ndir-1,2	,3		);
   &patdef('CE&',$ndir-1			);
   &patdef('CELOOP',$ndir-1			);

   &patdef('C'	,$ndir	,1	,2	,3	);
   &patdef('C&'	,$ndir				);
   &patdef('CLOOP',$ndir       			);
   &patdef('CC'	,$ndir	,'^C'	,'^C'	,'^C'	);

   &patdef('LIM'	,2	,'min'	,'max'	);
   &patdef('LLIM'	,2	,'lo'	,'hi'	);
   &patdef('L'		,2	,'min^D','max^D');
   &patdef('LL'		,2	,'lo^D'	,'hi^D');
   &patdef('LSUB'	,2	,'+'	,'-'	);
   &patdef('LADD'	,2	,'-'	,'+'	);
   &patdef('LT'		,2	,'>'	,'<'	);

   &patdef('LENTYPE'	,1	,10		);
   &patdef('LENNAME'	,1	,400	);

   &patdef('IFONED'	,$ndim==1		);
   &patdef('IFTWOD'	,$ndim==2		);
   &patdef('IFTHREED'	,$ndim==3		);
   &patdef('IFTHREEC'	,$ndir==3		);
   &patdef('NOONED'	,$ndim!=1		);
   &patdef('IFPHI'      ,$phi>0			);
   &patdef('IFZ'        ,$z>0			);
   &patdef('IFFCT'	,$if_fct  ,'      '	);
   &patdef('IFTVDLF'	,$if_tvdlf,'        '	);
   &patdef('IFTVD'	,$if_tvd  ,'      '	);
   &patdef('IFIMPL'	,$if_impl ,'       '	);
   &patdef('IFPOISSON'	,$if_poisson ,'        ');
   &patdef('IFCT'       ,$if_ct   ,'     '      );
   &patdef('IFGEN'	,$if_gencoord,'      '	);
   &patdef('NOGEN'	,!$if_gencoord,'      '	);
   &patdef('IFRES'	,$if_resist,'      '	);
   &patdef('NORES'	,!$if_resist,'      '	);
   &patdef('ANDIFRK'	,$if_rk   ,'        '   );
   &patdef('ANDIFCD'	,$if_cd   ,'        '	);
   &patdef('ANDIFMC'	,$if_mc   ,'        '	);
   &patdef('IFMPI'      ,$if_mpi  ,'     '      );
   &patdef('NOMPI'      ,!$if_mpi ,'     '      );
}
#============================================================================
sub patdef{
   # Put pattern definitions into global %sub and %nsub associative arrays
   local($pattern,$nsub,@substitute)=@_; local($isub);
   for($isub=0;$isub<=$#substitute;$isub++){
       $sub{$pattern,$isub+1}="$substitute[$isub]";
   }
   $nsub{$pattern}=$nsub;  
}
#============================================================================
sub block{
   # Find block opening by $ope and closed by matching $clo, and substitute
   local($input)=@_;
   local($rightpos,$leftpos,$line,$nline);
   $leftpos=index($_,$ope);
   $level=0;
   $rightpos=&rightmatch($leftpos,$ope,$clo,$clo);
   #TST print "Leftpos:$leftpos,Rightpos:$rightpos\n";
   # If no match, ie $rightpos>=length($_), read at most 100 more lines
   while($rightpos>=length($_) && $nline<100 && ($line=<$input>)){
       $_.=&unquote($line); $nline++;
       $rightpos=&rightmatch($rightpos-1,$ope,$clo,$clo);
       #TST print "Line:$_ Rightpos:$rightpos\n";
   }
   die"$ope without matching $clo within 100 lines" if $rightpos>=length($_);
   &substitute($leftpos,$rightpos);
}
#============================================================================
sub expression{
   # Find expression surrounding a bounded pattern and do substitution
   local($patpos,$rightpos,$leftpos);
   $patpos=index($_,'^');
   $leftpos = &leftmatch($patpos,'(',')',$lbound);
   $rightpos=&rightmatch($patpos,'(',')',$rbound);
   &substitute($leftpos,$rightpos);
}
#============================================================================
sub rightmatch{
   local($p,$lparen,$rparen,$bound)=@_; local($i,$c);
   for($i=$p+1;$i<length($_);$i++){
      $c=substr($_,$i,1); 
      last unless ($level!=0 || index($bound,$c)<0);
      $level++ if $c eq $lparen;
      $level-- if $c eq $rparen;
   }
   $i;
}
#============================================================================
sub leftmatch{
   local($p,$lparen,$rparen,$bound)=@_; local($i,$c,$level);
   for($i=$p-1;$i>0;$i--){
      $c=substr($_,$i,1); 
      last unless ($level!=0 || index($bound,$c)<0);
      $level++ if $c eq $rparen;
      $level-- if $c eq $lparen;
   }
   $i;
}
#============================================================================
sub substitute{
   #Substitute patterns of the same ID as of the first pattern
   local($leftpos,$rightpos)=@_;
   local($head,$lchr,$expr,$rchr,$tail);
   local($id,$nsub,$expa,$result);

   $head=substr($_,0,$leftpos);
   $lchr=substr($_,$leftpos,1);
   $expr=substr($_,$leftpos+1,$rightpos-$leftpos-1);
   $rchr=substr($_,$rightpos,1);
   $tail=substr($_,$rightpos+1);

   #Determine separator
   if($rchr eq ";")                 {$separator=";"}
   elsif($expr=~s/$bar([^$bar]*)$//){$separator="$spc$1$spc"}
   elsif($expr=~s/\\$//)            {$separator="$spc$brk$spc"}
   elsif($expr=~s/$sep$//)          {$separator="$&"}
   else                             {$separator=","};

   # Determine number of substitutions and ID
   if($expr=~/$pat($patid)($patchr*)/){$id=$1;$nsub=$nsub{"$1$2"};}
   else{$nsub=$nsubdefault};
   # Do the substitutions
   for($isub=1;$isub<=$nsub;$isub++){
       $expa=$expr;
       # If there is a unitvector pattern ^%N choose preceeding or subsequent
       # part depending on N==$isub or not.
       $expa=~s/$pat$uni$isub.*// || $expa=~s/.*$pat$uni.//;
       # Substitute patterns with $pat$id to their $isub-th substitute
       $expa=~s/$pat($id$patchr*)/$sub{$1,$isub}/g;
       $result.=$separator if $isub>1;
       $result.=$expa;
   }

   $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
   $rchr=';' if $rchr eq $clo && $separator eq ";";
   $rchr='' if $rchr eq $clo;
   $_=$head.$lchr.$result.$rchr.$tail;
}
#============================================================================
sub common{
    # Process lines of the form "COMMON, type:: var1(dim1),var2,..."
    # Store %common declarations with the name of the four first character
    # of the type and form the variable declaration type:: var1,var2,...
    /^ *COMMON, *(....)[^:]*:: */; $name=$1; $varlist=$';
    if($common{$name}){$common{$name}=~s/\n$/,$varlist/
    }else{             $common{$name}="COMMON /$name/ $varlist"}
    # Eliminate "COMMON," and the dimensions in parentheses
    s/^ *COMMON, *//; s/\([^)]*\)//g;
}
#============================================================================
sub shiftopt{
    # Process lines with special !SHIFT comments, optimize shifts.
    # !SHIFT and !SHIFT MORE is followed by a shift like 
    #     jxmin1=ixmin1+kr(idim,^D); ...
    # !SHIFT BEGIN and !SHIFT END enclose the lines where jx^L is used
    #     a(jxmin1:jxmax1,jxmin2:jxmax2)=b(ixmin1:ixmax1)...

    # Uses globals: $shiftdef, %shift0, %shiftby, 
    #               $shiftproc, $shiftregion, $shiftcount
    local($ss,$s0,$si,$sr,$sc,@cond,$shiftcond,@shifted,$shifted,
	  $idim,$indent);

    if($shiftdef){
        # Evaluate the shiftdefinition line: jxmin1 = ixmin1 - 3*kr...
        /(\w+) *= *(\w+) *([\+\-;]\d?)/
	    || die "Could not read shift in line $. !\n";
        $ss=$1; $s0=$2; $si=$3;
        # Determine shift from si="+2",";","-" or a similar string
        $si='' if $si=~/;/;
        $si.="1" if length($si)==1;
        # Store basis and shift into %shift0 and %shiftby without "min1"
        $ss=~s/min1//;$s0=~s/min1//;
        $shift0{$ss}=$s0;
        $shiftby{$ss}=$si;
        $shiftdef=0;
    }
    elsif(/!SHIFT *$/)   {$shiftdef=1; undef %shift0; undef %shiftby;}
    elsif(/!SHIFT MORE/) {$shiftdef=1}
    elsif(/!SHIFT BEGIN/){$shiftproc=1;}
    elsif(/!SHIFT END/){
        $indent=$`;
        $_='';
	die "Missing SHIFT BEGIN before $.\n" if $shiftregion eq '';
	#die "Missing SHIFT before $.\n" unless defined(%shift0);
	$shiftcount++;
        # Produce the Perl match string and the Fortran IF condition string
	@shifted=keys(%shift0);
	$shifted=join('|',@shifted);
	@shiftcond=@shifted;
	foreach $ss (@shiftcond){
	    for($idim=1;$idim<=$ndim;$idim++){
		$cond[$idim-1]=$ss."min$idim==".$shift0{$ss}."min$idim";
	    }
	    $ss=join('.and.',@cond);
	}
	$shiftcond=join(".and.\&\n$indent        ",@shiftcond);
	# OPTIMIZE SHIFT REGION
        for($idim=1;$idim<=$ndim;$idim++){
	    # CHECK IN FORTRAN IF THERE IS SHIFT IN DIRECTION IDIM
	    $sc=$shiftcond; 
	    foreach $ss (@shifted){
		$sc=~s/(${ss}min$idim==$shift0{$ss}min$idim)/$1$shiftby{$ss}/;
	    }
	    if($idim==1){$_.="${indent}IF     ($sc) THEN\n";}
	    else        {$_.="${indent}ELSE IF($sc) THEN\n";}
	    # REPLACE SHIFTED VARIABLES BY SHIFT0SHIFTBY IN THE SHIFT REGION
	    $sr=$shiftregion;
	    $sr=~s/\b($shifted)((min|max)$idim)\b/$shift0{$1}$2$shiftby{$1}/g;
            # IN OTHER DIRECTIONS THERE IS NO SHIFT SO USE JUST SHIFT0
            $sr=~s/\b($shifted)((min|max)[123])\b/$shift0{$1}$2/g;
	    $_.=$sr;
	}
        # OTHERWISE PRINT A WARNING AND SHIFT BY ORIGINAL STATEMENT
	$_.="${indent}ELSE\n$indent   IF(it==itmin)write(*,*)".
	    "'Warning in $filename: ',\&\n".
		"$indent      'SHIFT$shiftcount did not optimize!'\n";
	$_.="$shiftregion${indent}ENDIF\n$indent!SHIFT END\n";
        # CLEANUP
	undef $shiftregion;
	$shiftproc=0;
    }
    elsif($shiftproc){$shiftregion.="$_"; $_='';} # Collect lines
}
#============================================================================
sub printline{
    local($line,$comment); $sss="   ";

    # PRINT FORMATTED OUTPUT LINE BY LINE
    while(s/(.*)\n//){
       $line=$1;
       # PUT TRAILING COMMENTS INTO $comment
       if($line=~s/ *!.*//){
           $comment=&unquote($&);
           # If line is longer than $maxlen, try reducing length of comment
           if(length("$line$comment")>$maxlen){$comment=~s/ *! */ !/}
       }else{$comment=""};
       # Print full line
       $line=&unquote(&format90($line));
       # CORRECT REFERENCES TO INCLUDE FILES
       if($line=~/^ *include 'vac/){
          $line=~s/include 'vacdef\.\w+'/include 'vacdef.f'/;
          $line=~s/include 'vacpar\.\w+'/include 'vacpar.f'/;
          $line=~s/include 'vacusrpar\.\w+'/include 'vacusrpar.f'/;
       }
       print $line,"$comment\n";
    }
}
#===========================================================================
sub format90 {
   local($line)=@_;

   # If line is not too long return
   return($line) if length($line)<=$maxlen;

   # BREAK LONG LINES INTO CONTINUATION LINES AND/OR REDUCE INDENTATION
   local($bestlen,$goodlen,$len,$maxindent,$indent,$indentnow,$c,$answer);

   # Determine line indentation. If too much, reduce it to maximum.
   $maxindent=int(0.6*$maxlen);
   $line=~s/^( {0,$maxindent}) */$1/; $indent=$1;

   # We are happy if the length of the line is between $goodlen and $maxlen
   $goodlen=$maxlen-20;

   # Start breaking line
   while(length($line)>$maxlen){
       # Check for semicolon
       if(($len=rindex($line,';',$maxlen-1))>=0){
          $answer.=substr($line,0,$len)."\n";
          $line=substr($line,$len+1); $line=~s/^ */$indent/;
          next;
       }
       # Find best breakpoint after indentation but before $maxlen-2
       $line=~/^ */; 
       $bestlen=$indentnow=length($&);
       foreach $c ('+','-','//','/','**','*','.or.','.and.','>','<','==','=',
                   ',',' ','(',')'){
          if(($len=rindex($line,$c,$maxlen-length($c)))>=0){
             $len+=length($c) if $c =~ m!^,| |\)|//$!;
             $bestlen=$len if $len>$bestlen;
             last if $bestlen>$goodlen;
          }
       }
       if($bestlen>$indentnow){
          # Collect broken parts in $answer and break the continuation further
          $answer.=substr($line,0,$bestlen)."&\n";
          $line="$indent$sss".substr($line,$bestlen);
       }else{
          # No break was found. Remove indentation if there is any or die.
          $line=~s/^ +// || die "Couldn't break line:".&unquote($line) 
       }
   }
   $answer.=$line;
}
#===========================================================================
sub quotation{
   # Hide quoted strings by converting them from 0-127 to 128-255 ASCII.
   # Only matched quotation marks ( ' and " ) count.
   local($line)=@_; local($head,$q);

   # From left to right invert quoted text into $head. The rest is in $line.
   while($line=~/(['"])/){
      $head.=$`; $q=$1; $line=$';
      # Check for matching quotation mark. If found, convert quoted part.
      if($line=~/$q/){$q="$q$`$q"; $line=$'; $q=~tr/\x00-\x7f/\x80-\xff/;}
      $head.=$q;
   }
   # return result
   $head.$line;
   
}
#===========================================================================
sub unquote{
    # Uncover hidden quotations by deleting 8-th bit.
    local($line)=@_;
    $line=~tr/\x80-\xff/\x00-\x7f/;
    $line;
}
