#script to move along with the polymere, adding stars and trail.
#use for example:
#cat data.xyz | awk -f ~/Home/scripte/comovingstar3d.awk > temp.xyz ; vmd -e /Users/elgeti/Home/spp/spsemiflex/visu.vmd 

BEGIN{ #initialize stuff
   traill=100
    for (i=0;i<traill;i++){
	trailx[i]=25
	traily[i]=25
	trailz[i]=0
    }

stars=3
dist=40
xp=25
yp=25
zp=5.3
}
{ #begin main ---------------------------
if(NF==1){n=$1;line=0}  #number of particles in xyz file
line+=1
if(line==1) print n+8*stars*stars*stars+traill #n squared "stars"
if(line==2) {
    print $0
}
if(line==3) {
    #xp=$2;yp=$3;zp=$4  #position of first monomere
    #smother comotion:
    xp=xp+0.1*($2-xp)
    yp=yp+0.1*($3-yp)
    zp=zp+0.1*($4-zp)
#make trail:
    for (i=0;i<traill-1;i++){
	trailx[i]=trailx[i+1]
	traily[i]=traily[i+1]
	trailz[i]=trailz[i+1]
    }
    trailx[traill-1]=$2
    traily[traill-1]=$3
    trailz[traill-1]=$4

    for (i=0;i<traill;i++){
	print "S ",trailx[i]-xp,traily[i]-yp,trailz[i]-zp
    }    

    for (i = -stars; i < stars; i++) {  #place stars
	for (j = -stars; j < stars; j++) {  #place stars
	    for (k = -stars; k < stars; k++) {  #place stars
		print "H ", dist*(i+0.5)-xp+int(xp/dist)*dist, dist*(j+0.5)-yp+int(yp/dist)*dist,dist*(k+0.5)-zp+int(zp/dist)*dist
	    }}}
} 
if(line>2){
    print "C ",$2-xp,$3-yp,$4-zp
}

}