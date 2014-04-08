!=======================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!=======================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!=======================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!
!-----------------------------------------------------------------------
!
! SUBROUTINES IN THIS FILE WERE WRITTEN BY N. TEANBY UNIVERSITY OF LEEDS
! BASED ON TEANBY AND KENDALL (2003) (?)
!

c-----------------------------------------------------------------------
      subroutine zabs_max(a,n,np,amax)
c-----------------------------------------------------------------------
c
c      subroutine to return the maximum _absolute_ value of array a
c
c    in:
c      a(np)            real            array to be searched
c      n            int            number of data points
c      np            int            array dimension
c    out:
c      amax            real            maximum
c
c-----------------------------------------------------------------------
c      N. Teanby      6-8-02      Original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,i
      real a(np),amax

      amax=0.
      do 1 i=1,n
         amax = max(amax,abs(a(i)))
1      continue

      return
      end

c-----------------------------------------------------------------------
      subroutine zagglomerative_cluster(x,y,n,xc,yc,vxc,vyc,cluster)
c-----------------------------------------------------------------------
c
c      subroutine to perform agglomerative/hierarchical clustering on a
c      2D data set (x,y). Initially we have n clusters, at each step the
c      nearest two clusters are combined to give one less cluster. This
c      is continued until there is only one cluster, composed of the entire
c      dataset.
c
c      algorithm
c      ---------
c      do for no_cluster k = n - 1
c      -calc all inter-cluster distances (euclidian distance)
c      -find clusters A and B with minimum distance
c      -reassign points in cluster B to cluster A
c      -renumber clusters from 1-k
c      -calc mean and variance with each cluster
c      -calc the cluster criteria (maximise to find no. clusters)
c      repeat
c
c      variables
c    in:
c      x/y(npc)            real      x and y data to cluster
c      n                  int      number of points
c    out:
c      xc/yc(npc,npc)      real      cluster centres (=mean of datapoints in cluster)
c                              first index is cluster number (1-k)
c                              second index is number of clusters (=k=1-n)
c      vxc/vyc(npc,npc)      real      within cluster variance
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c
c    other:
c      k                  int      number of clusters
c      d(npc,npc)            real      distance between i and jth cluster
c
c      NB. matricies with two indicies are upper triangular. first index is
c      cluster number (1-k), second index is number of clusters (=k=1-n)
c
c-----------------------------------------------------------------------
c      n.teanby      20-08-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n
      real x(npc),y(npc)
      integer cluster(npc,npc),nc(npc)
      real xc(npc,npc),yc(npc,npc),vxc(npc,npc),vyc(npc,npc)
      real d(npc,npc),dmin
      integer i,k,imin,jmin

c  ** initially there are n clusters **
      k=n
c  ** assign the k clusters **
      do i=1,k
         cluster(i,k)=i
         xc(i,k)=x(i)
         yc(i,k)=y(i)
      enddo

c  ** reduce the number of clusters from n to 1 **
c  ** by grouping the nearest neigbours **
      do k=n-1,1,-1
c     ** calc the distance matrix between the k+1 clusters **
         call zdiss_euclid(xc,yc,k+1,k+1,npc,d)
c      ** find minimum distance **
         call zdiss_min(d,k+1,npc,imin,jmin,dmin)
c     ** reasign the datapoints in cluster imin to cluster jmin **
         do i=1,n
            if (cluster(i,k+1).eq.imin) then
               cluster(i,k)=jmin
            else
               cluster(i,k)=cluster(i,k+1)
            endif
         enddo
c     ** renumber clusters from 1-k (i.e. remove gaps in the cluster nos) **
         call zclust_renumber(cluster,n,npc,k)
c     ** find the average cluster positions **
         call zcluster_loc(x,y,cluster,k,n,npc,xc,yc,vxc,vyc,nc)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine zclust_renumber(cluster,n,np,k)
c-----------------------------------------------------------------------
c
c      renumber an array of k +ve integers so that they number from 1-k
c
c      uses a collapsing algorithm where the excess differences between
c      elements are incrementally removed
c
c      cluster(np,np)      int      cluster numbers
c      n                  int      number of data
c      np                  int      array dimension
c      k                  int      number of distinct clusters
c
c-----------------------------------------------------------------------
c      n.teanby      9-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,jump,imax
      integer cluster(np,np)
      integer i,j,k

c  ** find maximum cluster number **
      imax=cluster(1,k)
      do i=1,n
         if (cluster(i,k).gt.imax) then
            imax=cluster(i,k)
         endif
      enddo

c  ** renumber the clusters from 1 to k **
c  ** for every cluster **
      do i=1,k
         jump=imax
c      ** for all the data points **
         do j=1,n
            if (cluster(j,k).eq.i) then
c            ** if cluster number exists move to next cluster number **
               goto 1
            else
c            ** else find the difference between cluster no. i and the nearest
c               old cluster number **
               if (cluster(j,k).gt.i) then
                  jump=min0(jump,cluster(j,k)-i)
               endif
            endif
         enddo
c      ** shift all clusters numbered over i by -jump **
         do j=1,n
            if (cluster(j,k).gt.i) then
               cluster(j,k)=cluster(j,k)-jump
            endif
         enddo
1         continue
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine zdiss_euclid(x,y,k,n,np,d)
c-----------------------------------------------------------------------
c
c      calculate the Euclidian distances (d=sqrt(dx**2+dy**2)) between a
c      set of n points defined by the kth column of x and y
c
c      x/y(np,np)            int      x/y coordinates of points
c      k                  int      column to use in calculation
c      n                  int      number of data
c      np                  int      array dimension
c      d(np,np)            int      distance between i and j th point
c
c      NB. d is only defined for elements i=2,n and j=1,i-1 (ie lower
c      triangular matrix). diagonal is undefined as it reperesents the
c      distance between the same point (=0)
c
c-----------------------------------------------------------------------
c      n.teanby      9-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,k
      real x(np,np),y(np,np),d(np,np)
      integer i,j

      do i=2,n
         do j=1,i-1
            d(i,j)=sqrt((x(i,k)-x(j,k))**2+(y(i,k)-y(j,k))**2)
         enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine zdiss_min(d,n,np,imin,jmin,dmin)
c-----------------------------------------------------------------------
c
c      search the dissimalarity matrix d for the lowest value
c
c      d(np,np)            real      dissimalarity matrix
c      n                  int      number of data
c      np                  int      array dimension
c      imin                  int      i location of minimum dissimalarity
c      jmin                  int      j location of minimum dissimalarity
c      dmin                  real      minimum dissimalarity =d(imin,jmin)
c
c      NB. d is only defined for elements i=2,n and j=1,i-1 (ie lower
c      triangular matrix)
c      jmin < imin (because j < i)
c-----------------------------------------------------------------------
c      n.teanby      9-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,imin,jmin
      real d(np,np),dmin
      integer i,j

c  ** set initial trial minimum **
      dmin=d(2,1)
      imin=2
      jmin=1
c  ** find minimum distance **
      do i=2,n
         do j=1,i-1
            if (d(i,j).lt.dmin) then
               dmin=d(i,j)
               imin=i
               jmin=j
               endif
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine zcluster_loc(x,y,cluster,k,n,np,xc,yc,vxc,vyc,nc)
c-----------------------------------------------------------------------
c
c      find the average cluster positions, variances, and number in each
c      cluster
c
c      x/y(np)            real      x/y possitions of data points
c      cluster(np,np)      int      k th column gives cluster number of i th data point
c      k                  int      number of cluster (column to use in calculation)
c      n                  int      number of data
c      np                  int      array dimension
c      xc/yc(np,np)      real      mean location of i th cluster when there
c                              are a total ofk clusters
c      vxc/vyc(np,np)      real      variance in location of i th cluster
c      nc(np)            int      number of data in i th cluster
c
c      NB. (v)x/yc is only defined for elements i=2,n and j=1,i-1 (ie lower
c      triangular matrix)
c-----------------------------------------------------------------------
c      n.teanby      9-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,k
      real x(np),y(np),xc(np,np),yc(np,np),vxc(np,np),vyc(np,np)
      integer cluster(np,np),nc(np)
      real xsum,ysum
      integer i,j

c  ** calc the mean positions **
      do j=1,k
         nc(j)=0
         xsum=0.
         ysum=0.
         do i=1,n
            if (cluster(i,k).eq.j) then
                  xsum = xsum + x(i)
                  ysum = ysum + y(i)
                  nc(j) = nc(j) + 1
               endif
         enddo
         xc(j,k)=xsum/real(nc(j))
         yc(j,k)=ysum/real(nc(j))
      enddo

c  ** calc the within cluster variance for each cluster **
      do j=1,k
         nc(j)=0
         xsum=0.
         ysum=0.
         do i=1,n
            if (cluster(i,k).eq.j) then
               xsum = xsum + (x(i)-xc(j,k))**2
               ysum = ysum + (y(i)-yc(j,k))**2
                  nc(j) = nc(j) + 1
               endif
         enddo
         vxc(j,k)=xsum/real(nc(j))
         vyc(j,k)=ysum/real(nc(j))
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine zwrite_clustxy(lu,file_clustxy,n,np,
     >      wbeg,wend,fast,dfast,tlag,dtlag)
c-----------------------------------------------------------------------
c      write out measurements to file
c-----------------------------------------------------------------------
c      n.teanby      9-9-02      original code
c-----------------------------------------------------------------------

      implicit none
      integer lu,n,np,i
      real wbeg(np),wend(np),fast(np),dfast(np),tlag(np),dtlag(np)
      character*50 file_clustxy

c  ** write file **
      open(lu,file=file_clustxy,status='unknown')
      do i=1,n
      write(lu,100) i,wbeg(i),wend(i),fast(i),dfast(i),tlag(i),dtlag(i)
      enddo
      close(lu)
100   format(i6,2f12.4,f8.3,f7.3,2f10.6)
      return
      end

c-----------------------------------------------------------------------
      subroutine zwrite_clusters(lu,file_clusters,n,np,
     >      xc0,yc0,vxc0,vyc0)
c-----------------------------------------------------------------------
c      write out clusters to file
c-----------------------------------------------------------------------
c      n.teanby      9-9-02      original code
c-----------------------------------------------------------------------

      implicit none
      integer lu,n,np,i
      real xc0(np),yc0(np),vxc0(np),vyc0(np)
      character*50 file_clusters

c  ** write file **
      open(lu,file=file_clusters,status='unknown')
      do i=1,n
         write(lu,*) xc0(i),yc0(i),sqrt(vxc0(i)),sqrt(vyc0(i))
      enddo
      close(lu)
      return
      end

c-----------------------------------------------------------------------
      subroutine zass_ini(inifile,lu,ext1,ext2,nwbeg,nwend,
     >dt_beg,dt_end,dtlag_max,dfast_max,t_off_beg,t_off_end,
     >tlag_scale,fast_scale,max_no_clusters,nmin,
     >OPT_verbose,OPT_outfiles)
c-----------------------------------------------------------------------
c
c      read in fast_scale and tlag_scale from the ass.ini file
c
c      variables
c      ---------
c    in:
c      inifile      char*50      name of parameters file
c      lu            int            logical unit to open files on
c    out:
c      ext1/2      char*50      extention of the two components to do analysis on
c      nwbeg            int            number of start positions for S-wave window
c      nwend            int            number of end positions for S-wave window
c      dt_beg      real            increment for start of window (in seconds)
c      dt_end      real            increment for end of window (in seconds)
c      t_off_beg      real            maximum window beginning (relative to spick)
c      t_off_end      real            minimum window end (relative to spick)
c      dtlag_max      real            max allowable error in lag time for inclusion
c                              in clustering
c      dfast_max      real            " fast direction "
c      tlag_scale      real            range of tlag scale in seconds
c      fast_scale      real            range of fast direction scale in degrees
c      max_no_clusters      int      max. number of clusters
c      nmin            int            minimum number of points in an acceptable cluster
c      OPT_verbose      logical      true for verbose output
c      OPT_outfiles logical      true if write outfiles for gmt plots
c
c-----------------------------------------------------------------------
c      n.teanby      12-5-03      original code
c-----------------------------------------------------------------------

      implicit none
      integer nwbeg,nwend,lu
      integer nmin,max_no_clusters
      real dtlag_max,dfast_max,tlag_scale,fast_scale,dt_beg,dt_end
      real t_off_beg,t_off_end
      logical OPT_verbose,OPT_outfiles
      character*50 inifile,ext1,ext2

      open(lu,file=inifile,status='old')
      read(lu,*) ext1
      read(lu,*) ext2
      read(lu,*) nwbeg
      read(lu,*) nwend
      read(lu,*) dt_beg
      read(lu,*) dt_end
      read(lu,*) dtlag_max
      read(lu,*) dfast_max
      read(lu,*) t_off_beg
      read(lu,*) t_off_end
      read(lu,*) tlag_scale
      read(lu,*) fast_scale
      read(lu,*) max_no_clusters
      read(lu,*) nmin
      read(lu,*) OPT_verbose
      read(lu,*) OPT_outfiles
      close(lu)

      return
      end
c-----------------------------------------------------------------------
      subroutine zbreadsac(file,lu,np,h1,h2,h3a,h3b,h4,d1,npts)
c-----------------------------------------------------------------------
c
c      read in a sac binary data file
c        read in header info (in 4 parts)
c        read in data part1
c        (may need to add a data part2 if unevenly spaced data is ever used)
c
c      variables:
c input: file      char*50            file to read in
c         lu            int                  logical unit to read file to
c         np            int                  dimension of d1 array
c output:h1            real(14,5)            header part1
c         h2            int(8,5)            header part2
c         h3a      char*8            header part3a
c         h3b      char*16            header part3b
c         h4            char*8 (7,3)      header part4
c         d1            real(np)            data points
c         npts      int                  number of data points (read from h2(2,5))
c
c-----------------------------------------------------------------------
c  modifications:
c      17-09-02      N. Teanby      Original code
c-----------------------------------------------------------------------

      implicit none

      integer lu,np,npts
      real h1(14,5)
      integer h2(8,5)
      character h3a*8, h3b*16
      character*8 h4(7,3)
      real d1(np)
      character*4 char_tmp1,char_tmp2,char_tmp3,char_tmp4
      character*50 file
      integer i,j

c  ** open file as direct access **
      open(unit=lu,file=file,form='unformatted',access='direct',recl=4,
     >status='old')

c  ** h1 **
      do i=1,14
         do j=1,5
            read(lu,rec=5*(i-1)+j) h1(i,j)
         enddo
      enddo
c  ** h2 **
      do i=1,8
         do j=1,5
            read(lu,rec=70+5*(i-1)+j) h2(i,j)
         enddo
      enddo
c  ** h3 **
       read(lu,rec=111) char_tmp1
       read(lu,rec=112) char_tmp2
      h3a=char_tmp1//char_tmp2
       read(lu,rec=113) char_tmp1
       read(lu,rec=114) char_tmp2
       read(lu,rec=115) char_tmp3
       read(lu,rec=116) char_tmp4
      h3b=char_tmp1//char_tmp2//char_tmp3//char_tmp4
c  ** h4 **
      do i=1,7
         do j=1,3
            read(lu,rec=116+6*(i-1)+2*j-1) char_tmp1
            read(lu,rec=116+6*(i-1)+2*j) char_tmp2
            h4(i,j)=char_tmp1//char_tmp2
         enddo
      enddo
c  ** d1 **
c  ** read npts from header **
      npts=h2(2,5)
c  ** check npts does not exceed array dimension **
      if (npts.gt.np) then
         print*,'WARNING:zsacread'
         print*,'WARNING:  npts.gt.np, data missed out'
      endif
c  ** read in data **
      do i=1,npts
         read(lu,rec=158+i) d1(i)
      enddo

      close(lu)
      end
c-----------------------------------------------------------------------
      subroutine zc_73dudhar(x,y,n,xc,yc,xmin,ymin,
     >cluster,max_no_clusters,c,kopt)
c-----------------------------------------------------------------------
c
c      calc the clustering criteria of Duda and Hart 1973
c
c      consider 2 groups which are combined to form a single group in the
c      next iteration.
c
c      je2 is the sum of square errors within the two clusters
c      je1 is the sum of square errors within the combined cluster
c
c      this rule is for continuous data and we have discrete data so
c      the minimum errors are set to xmin/ymin, corresponding to the grid
c      spacing. this also avoids division by zero issues
c
c      calc je2/je1
c      if je2/je1 > c_critical then stop the clustering here
c
c      Milligan and Cooper 1985 found c_critical = 3.20 gave good results
c
c      variables
c      ---------
c    in:
c      x/y(npc)            real      data that has been clustered
c      n                  real      number of data points
c      xc/yc(npc,npc)      real      cluster centres (=mean of datapoints in cluster)
c                              first index is cluster number (1-k)
c                              second index is number of clusters (=k=1-n)
c      x/ymin            real      grid spacing of x/y data
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c      max_no_cluster      int      max number of clusters
c    out:
c      c(npc)            real      clustering criteria to determine optimum
c                              number of clusters
c      kopt                  int      optimum number of clusters
c
c-----------------------------------------------------------------------
c      n.teanby      28-08-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,kopt,max_no_clusters
      real x(npc),y(npc),xc(npc,npc),yc(npc,npc),c(npc)
      integer cluster(npc,npc)
      real je1,je2,c_critical,xmin,ymin
      integer i,j,k,cluster1,cluster2,cluster12,nc
      logical first_pass

      kopt=0
      c(n)=0.
      c_critical=3.20

c  ** for all numbers of clusters **
      do k=1,n-1
c      ** for each cluster **
         do j=1,k
c         ** for each data point **
c         ** search each data point to find a single cluster which is two seperate clusters when the number of clusters is increased by one **
            first_pass=.true.
            do i=1,n
               if (cluster(i,k).eq.j) then
                  if (first_pass) then
                     first_pass=.false.
                     cluster1 = cluster(i,k+1)
                     cluster12= j
                  else if (cluster(i,k+1).ne.cluster1) then
                     cluster2=cluster(i,k+1)
                     goto 99
                  endif
               endif
            enddo
         enddo
99         continue
c  ** calc Je1 **
c  ** cluster12 = number of the combined cluster when there are k clusters **
         je1=0.
         nc=0
         do i=1,n
            if (cluster(i,k).eq.cluster12) then
                  je1 = je1 + max(xmin,x(i)-xc(cluster12,k))**2 +
     >                     max(ymin,y(i)-yc(cluster12,k))**2
               nc=nc+1
            endif
         enddo
c  ** calc Je2 **
c  ** cluster1 and cluster2 = number of the seperate clusters when there
c      are k+1 clusters **
         je2=0.
         do i=1,n
            if (cluster(i,k+1).eq.cluster1) then
                  je2 = je2 + max(xmin,x(i)-xc(cluster1,k+1))**2 +
     >                     max(ymin,y(i)-yc(cluster1,k+1))**2
            endif
            if (cluster(i,k+1).eq.cluster2) then
                  je2 = je2 + max(xmin,x(i)-xc(cluster2,k+1))**2 +
     >                     max(ymin,y(i)-yc(cluster2,k+1))**2
            endif
         enddo
c      ** this is the special case for 2 parameters **
         c(k)=(0.681690113 - je2/je1)*sqrt(real(nc)/0.18943053)
      enddo

c  ** search for the optimum number of clusters **
c  ** if c(k) exceeds c_critical then the cluster should be subdivided,
c      giving kopt = k + 1 as the optimum number of clusters **
      do k=1,max_no_clusters
         if (c(k).gt.c_critical) then
            kopt = max(kopt,k+1)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zc_74calhar(x,y,n,xc,yc,xmin,ymin,
     >cluster,max_no_clusters,c,kopt)
c-----------------------------------------------------------------------
c
c      calc the clustering criteria of Calinski and Harabasz 1974
c
c             trace(B)/(k-1)
c      c(k) = --------------
c             trace(W)/(n-k)
c
c      W = within-cluster covariance matrix
c      B = between-cluster covariance matrix
c      k = number of clusters
c      n = number of datapoints
c
c      maximise c(k) to find the optimum number of clusters
c
c      this rule is for continuous data and we have discrete data so
c      the minimum errors are set to xmin/ymin, corresponding to the grid
c      spacing. this also avoids division by zero issues
c
c      variables
c      ---------
c    in:
c      x/y(npc)            real      data that has been clustered
c      n                  real      number of data points
c      xc/yc(npc,npc)      real      cluster centres (=mean of datapoints in cluster)
c                              first index is cluster number (1-k)
c                              second index is number of clusters (=k=1-n)
c      x/ymin            real      grid spacing of x/y data
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c      max_no_cluster      int      max number of clusters
c    out:
c      c(npc)            real      clustering criteria to determine optimum
c                              number of clusters
c      kopt                  int      optimum number of clusters
c
c-----------------------------------------------------------------------
c      n.teanby      28-08-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,kopt,max_no_clusters
      real x(npc),y(npc),xc(npc,npc),yc(npc,npc),c(npc)
      integer cluster(npc,npc)
      real xbar,ybar,traceb,tracew,c_max,xmin,ymin
      integer i,j,k,nc

c  ** calc x/ybar **
      xbar=0.
      ybar=0.
      do i=1,n
         xbar=xbar+x(i)
         ybar=ybar+y(i)
      enddo
      xbar=xbar/real(n)
      ybar=ybar/real(n)

c  ** calc c for each number of clusters **
      c(1)=0.
      do k=2,n
         tracew=0.
         traceb=0.
         do j=1,k
            nc=0
            do i=1,n
               if (cluster(i,k).eq.j) then
                  tracew=tracew+(max(xmin,x(i)-xc(j,k)))**2
     >                        + (max(ymin,y(i)-yc(j,k)))**2
                  nc=nc+1
               endif
            enddo
            traceb=traceb+nc*((xc(j,k)-xbar)**2+(yc(j,k)-ybar)**2)
         enddo
         c(k)=traceb*real(n-k)/(real(k-1)*tracew)
c         print*,'74calhar',k,traceb,tracew,c(k)
      enddo

c  ** find optimum number of clusters **
      c_max=c(1)
      kopt=1
      do k=2,max_no_clusters
         if (c(k).gt.c_max) then
            c_max=c(k)
            kopt = k
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zcluster(x0,dx0,y0,dy0,n,
     >xscale,yscale,xmin0,ymin0,max_no_clusters,
     >xc0,yc0,vxc0,vyc0,cluster,k)
c-----------------------------------------------------------------------
c
c      subroutine to perform cluster analysis on a dataset comprised of
c      tlag and fast direction from the different trial windows. Cluster
c      the data and find the optimum number of clusters.
c
c    in:
c      x0/y0(npc)             real      tlag and fast values
c      dx0/dy0(npc)      real      standard deviation of tlag and fast measurements
c      n                  int      number of data points
c      x/yscale            real      scale/standardisation factors for x0/y0 data
c      x/ymin0            real      grid spacing of x/y data
c      max_no_cluster      int      max number of clusters
c    out:
c      k                  int      optimum number of clusters = max(k1,k2)
c      xc0/yc0(npc)      real      cluster means
c      vxc0/vyc0(npc)      real      within cluster variance
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c    other:
c      x/y/dx/dy(npc)      real      scaled data and standard deviation
c      xc/yc(npc,npc)      real      scaled cluster means of ith cluster when k=j
c                               [indecies are (i,j)]
c      vxc/vyc(npc,npc)      real      scaled cluster variance of ith cluster when k=j
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c      c(npc)            real      clustering criteria to determine optimum
c                              number of clusters
c      k1                  int      optimnum cluster no. from Calinski and Harabasz
c                              1974 method
c      k2                  int      optimnum cluster no. from Duda and Hart
c                              1973 method
c-----------------------------------------------------------------------
c      n.teanby      15-8-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,j,k,k1,k2
      real x0(npc),y0(npc),dx0(npc),dy0(npc)
      real xc(npc,npc),yc(npc,npc),vxc(npc,npc),vyc(npc,npc)
      real xc0(npc),yc0(npc),vxc0(npc),vyc0(npc)
      real x(npc),y(npc),dx(npc),dy(npc)
      real xscale,yscale,xmin,ymin,xmin0,ymin0,c(npc)
      integer cluster(npc,npc)
      integer max_no_clusters

c  ** scale the data **
      call zmultiply(x0,n,npc,1./xscale,x)
      call zmultiply(y0,n,npc,1./yscale,y)
      call zmultiply(dx0,n,npc,1./xscale,dx)
      call zmultiply(dy0,n,npc,1./yscale,dy)
      xmin=xmin0/xscale
      ymin=ymin0/yscale

c  ** calc the cluster hierachy of the scaled data **
      call zagglomerative_cluster(x,y,n,xc,yc,vxc,vyc,cluster)

c  ** calc clustering criteria **
c  ** calc k using Calinski and Harabasz (1974) method **
      call zc_74calhar(x,y,n,xc,yc,xmin,ymin,
     >cluster,max_no_clusters,c,k1)
c  ** calc k using Duda and Hart (1973) method **
      call zc_73dudhar(x,y,n,xc,yc,xmin,ymin,
     >cluster,max_no_clusters,c,k2)

c  ** set number of clusters to max of k1 and k2 **
      k=max(k1,k2)

c  ** unscale the cluster possitions and variances **
      do j=1,k
         xc0(j) = xscale * xc(j,k)
         yc0(j) = yscale * yc(j,k)
         vxc0(j) = xscale**2 * vxc(j,k)
         vyc0(j) = yscale**2 * vyc(j,k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zcovariance(x,y,n,np,cov)
c-----------------------------------------------------------------------
c
c      calculate the 2x2 covariance matrix of the vectors x and y
c      variables
c      x(np)                        real      time series
c      y(np)                        real      time series
c      n                        int      number of points
c      np                        int      array dimension
c      cov(2,2)                  real      covariance matrix
c
c-----------------------------------------------------------------------
c      N. Teanby      16-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i
      real x(np),y(np),cov(2,2)

c  ** calc cross correlation matrix **
      cov(1,1) = 0.
      cov(2,2) = 0.
      cov(1,2) = 0.
      cov(2,1) = 0.
      do 1 i=1,n
         cov(1,1) = cov(1,1) + x(i)**2
         cov(2,2) = cov(2,2) + y(i)**2
         cov(1,2) = cov(1,2) + x(i)*y(i)
1      continue
      cov(2,1) = cov(1,2)

c  ** normalise **
ccc      cov(1,1) = cov(1,1)/real(n)
ccc      cov(1,2) = cov(1,2)/real(n)
ccc      cov(2,2) = cov(2,2)/real(n)

      return
      end

c-----------------------------------------------------------------------
      subroutine zdetrend(y,n,np,ydetrend)
c-----------------------------------------------------------------------
c
c      subroutine to remove trend from a dataset using a least squares
c      line of best fit (see squires p38 for theory)
c
c      variables
c    in:
c      y(np)            real            series to detrend
c      n            int            number of data points
c      np            int            array dimension
c    out:
c      ydetrend(np)real            detrended dataset
c    local:
c      m            real            grad of line of best fit
c      c            real            intercept of line of best fit
c
c-----------------------------------------------------------------------
c      N. Teanby      5-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i
      real y(np),ydetrend(np)
      real xmean,ymean,sumy,sum1,sum2,trend,m,c,TINY
      parameter (TINY=1.e-10)

      sum1=0.
      sum2=0.
      sumy=0.

c  ** calc mean of x, which goes from 1 to n, so mean is **
      xmean=real(n+1)/2.
c  ** calc the mean of y **
      do 1 i=1,n
         sumy = sumy + y(i)
1      continue
      ymean=sumy/real(n)

c  ** calc grad of line according to formulat in squires **
      do 2 i=1,n
         sum1 = sum1 + (real(i)-xmean)*y(i)
         sum2 = sum2 + (real(i)-xmean)**2
2      continue

c  ** calc grad and intercept **
      if (sum1.lt.TINY) then
         m = 0.
         c = ymean
      else
        m = sum1/sum2
        c = ymean - m*xmean
      endif

c  ** dtrend the data **
      do 3 i=1,n
         trend = m*real(i) + c
         ydetrend(i) = y(i) - trend
3      continue

      return
      end
c-----------------------------------------------------------------------
      subroutine zeigen2x2(matrix,lambda1,lambda2,vec1,vec2)
c-----------------------------------------------------------------------
c
c      calculate eigenvalues and eigenvectors of a 2x2 matrix
c      variables
c    in:
c      matrix(2,2)                  real      matrix to find eigenvalues of
c    out:
c      lambda1/2                  real      eignevalues (lambda1 is largest)
c      vec1/2                  real      eignevectors (vec1 corresponds to lambda1)
c
c      see Boas p413 for maths
c
c-----------------------------------------------------------------------
c      N. Teanby      16-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      real matrix(2,2),lambda1,lambda2,a,b,c,temp,norm,vec1(2),vec2(2)

c  ** EIGENVALUES **
c  ** eigenvalue are the solution of a quadratic eqn with c.f.s **
      a = 1.
      b = - matrix(1,1) - matrix(2,2)
      c = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)

      lambda1 = 0.5*( -b + sqrt(b**2 - 4*a*c))
      lambda2 = 0.5*( -b - sqrt(b**2 - 4*a*c))

c  ** order the eigenvalues so that lambda2 is the smallest **
      if (lambda2.gt.lambda1) then
         temp    = lambda1
         lambda1 = lambda2
         lambda2 = temp
      endif

C  ** EIGENVECTORS (normalised) **
      vec1(1)=1.
      vec1(2)=(lambda1-matrix(1,1))/matrix(1,2)
      norm = sqrt(vec1(1)**2 + vec1(2)**2)
      vec1(1)=vec1(1)/norm
      vec1(2)=vec1(2)/norm

      vec2(1)=1.
      vec2(2)=(lambda2-matrix(1,1))/matrix(1,2)
      norm = sqrt(vec2(1)**2 + vec2(2)**2)
      vec2(1)=vec2(1)/norm
      vec2(2)=vec2(2)/norm

      return
      end
c-----------------------------------------------------------------------
      subroutine zerror95(error,ndf,lambda2_min,ierror,jerror)
c-----------------------------------------------------------------------
c
c      subroutine to calculate normalise contours so that 1=95% confidence
c      contour. Also calc the errorbars on tlag and fast direction
c
c      calculated confidence interval according to appendix in
c      Silver and Chan 1991 using the formula:
c
c      lambda2    <=    1   +   k   f_(k,ndf-k) (1-alpha)
c      -------               -------
c      lambda2_min           ndf - k
c
c      where:
c      k = number of parameters = 2
c      alpha = confindence level = 0.05 for 95% confidence
c      f is the inverse of the F distribution
c      ndf = number of degrees of freedom
c      lambda2 = second eigen value of covaraiance matrix (= elements of error)
c      lambda2_min = minimum lambda2
c
c      calc errorbars on tlag and fast by finding the bounding rectangle of the
c      95% confidence contour.
c
c      variables
c    in:
c      error(np1,np2)            real      array of lambda2 (altered by routine)
c      np1/2                        int      array dimension
c      ndf                        int      number of degrees of freedom
c      lambda2_min                  real      minimum lambda2 (corresponds to solution)
c    out:
c      error(np1,np2)            real      returned normalised array of lambda2
c      ierror                  real      extent (half width) of contour in i dirn
c      jerror                  real      extent (half width) of contour in j dirn
c    other:
c      lambda2_max                  real      lambda2 corresponding to 95% confidence
c
c-----------------------------------------------------------------------
c      N. Teanby      1-8-02      Original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------

      implicit none
      integer ndf
      integer i,j,jmin,jmax,irange,jrange
      integer line(npc)
      real ierror,jerror
      integer k,k1,irange_min,irange_max,istart,line_test(npc)
      real error(np1,np2int),lambda2_min,fftable,lambda_max
      external fftable

c  ** check that npc is beg enough **
      if (npc.lt.np1) then
         pause 'ERROR: zerror95: npc < np1'
      endif

c  ** calc value of lambda at 95% confidence limit from tabulated values **
      if (ndf.ge.3) then
         lambda_max = lambda2_min*( 1. + 2.*fftable(ndf-2)/real(ndf-2))
      else
         print*,'WARNING: zerror95: ndf <=2, set to 3)'
         lambda_max = lambda2_min*( 1. + 2.*fftable(1))
      endif

c  ** normalise errors by lambda_max, so that and error of 1 = 95% confidence **
      do i=1,np1
         do j=1,np2int
            error(i,j)=error(i,j)/lambda_max
         enddo
      enddo

c  ** find i/j error (half width of 95% confidence contour) **
c  ** find min and max j, simply search the array **
      jmin=np2int
      jmax=1
      do i=1,np1
         do j=1,np2int
            if (error(i,j).le.1.0) then
               jmin = min0(jmin,j)
               jmax = max0(jmax,j)
            endif
         enddo
      enddo
      jrange=jmax-jmin
c  ** finding min max i is more difficult because of cyclicity of angles **
c  ** sweep a line over all j, set point on line equal to 1 if it falls within the 95% convidence contour for any j. The height of the bounding rectangle is defined by the shortest line which includes all points with a value of line(i)=1. This line is found by searching all line lengths from the minimum = sum_i linr(i) to maximum = np1**
      do i=1,np1
         line(i)=0
      enddo
      do j=1,np2int
         do i=1,np1
            if (error(i,j).le.1.0) then
               line(i)=1
            endif
         enddo
      enddo
c  ** min line length **
      irange_min = 0
      do i=1,np1
         irange_min = irange_min + line(i)
      enddo
c  ** max line length **
      irange_max = np1
c  ** search all line length and starting points to find irange **
      do i = irange_min , irange_max
        do istart = 1,np1
          do k = 1,np1
            line_test(k)=0
          enddo
          do k = istart,istart+i
            if (k.gt.np1) then
               k1 = k - np1
            else
               k1 = k
            endif
            line_test(k1) = 1
          enddo
          do k = 1,np1
            if ((line(k).eq.1).and.(line_test(k).ne.1)) then
              goto 1
            endif
          enddo
          irange = i
          goto 11
1          continue
        enddo
      enddo
11      continue

c  ** one standard deviation = 0.5*95% error half width
c      (so x full width of 95% contour by 0.25 to get 1s.d.)**
      ierror = 0.25*real(irange)
      jerror = 0.25*real(jrange)

      return
      end
c-----------------------------------------------------------------------
      subroutine zerror_interp(error,error_int)
c-----------------------------------------------------------------------
c
c      interpolate the error surface in the tlag direction.
c      the interpolation is done one row at a time
c
c      variables
c    in:
c      error(np1,np2)      real            error surface (i.e. lambda2)
c      np1/2                  int            array dimensions
c      np2int            int            np2 after interpolation
c    out:
c      error_int(np1,np2int)      real      interpolated error surface
c    local:
c      np                  int            array dimension
c                              (read from SIZE_np.h at compile time)
c
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      real error(np1,np2),error_int(np1,np2int)
      real error_row(np),error_row_int(np)
      integer f,i,j,n_int

c  ** check that np is big enough **
      if (np.lt.np2int) then
         pause 'ERROR: zerror_interp: np not big enough'
      endif

c  ** calc the interpolation factor based on np2 and np2int **
      f = (np2int-1)/(np2-1)
c  ** f needs to be an integer for zsplint to work **
      if (mod(np2int-1,np2-1).ne.0) then
         pause 'ERROR: zerror_interp: f not a whole number'
      endif

c  ** interpolate error surface in tlag direction **
c  ** do interpolation one row at a time **
      do 1 i=1,np1
c      ** copy ech row to a dummy array **
         do 11 j=1,np2
            error_row(j)=error(i,j)
11         continue
c      ** interpolate the row data **
         call zsplint(error_row,np2,f,error_row_int,n_int)
c     ** check that n_int = np2int
c         (this should not be possible but check anyway)**
         if (np2int.ne.n_int) then
            pause 'ERROR: zerror_interp: np2int.ne.n_int'
         endif
c      ** copy interpolated row to output array **
         do 22 j=1,n_int
            error_int(i,j)=error_row_int(j)
22         continue
1      continue

      return
      end
c-----------------------------------------------------------------------
      subroutine zerror_min(error,np1,np2,ifast,itlag,lambda2_min)
c-----------------------------------------------------------------------
c
c      grid search for the minimum value of lambda2 on the interpolated
c      error surface
c
c      variables
c    in:
c      error(np1,np2)      real      error surface (i.e. lambda2)
c      np1/2                  int      array dimensions
c    out:
c      lambda2_min            real      minimum value of lambda2
c      itlag                  int      index of lag corresponding to lambda2_min
c      ifast                  int      index of fast dirn corresponding to lambda2_min
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer np1,np2
      real error(np1,np2),lambda2_min
      integer i,j,itlag,ifast

c  ** find the minimum lambda2 position **
      lambda2_min=error(1,1)
      itlag=1
      ifast=1
      do 1 i=1,np1
        do 2 j = 1,np2
           if (error(i,j).lt.lambda2_min) then
              lambda2_min = error(i,j)
              itlag = j
              ifast= i
           endif
2         continue
1      continue

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE zfft(y,nn,np,isign,yfft)
c-----------------------------------------------------------------------
c
c      routine to calc fft of of a real time series
c
c      INPUT VARIABLES
c      y(np)            real            time seris
c      nn            int            number of data points
c      np            int            size of array
c      isign            int            set =  1 for fft
c                              set = -1 for inverse fft
c      OUTPUT VARIABLES
c      yfft(2*np)      real            fft of y
c            -odd and even indices are real and imag parts, respectively
c            -contains nn pairs of points
c
c      NOTES: nn must be a power of 2 (not checked for in routine)
c
c-----------------------------------------------------------------------
c      n. teanby      5-7-02      modified from Num Rec four1.f
c-----------------------------------------------------------------------
      INTEGER isign,nn,np
      REAL y(np),yfft(2*np)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp

c  -- this bit added --
c  ** copy y to yfft, because NR routine replaces y with it's fft **
c  ** yfft: odd and even indices are real and imag parts, respectively **
      do 100 i=1,nn
c      ** real part **
         yfft(2*i-1)    = y(i)
c      ** imaginary part (=0 because real data) **
         yfft(2*i)  = 0.0
100   continue

c  -- this bit unchanged from NR routine --
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=yfft(j)
          tempi=yfft(j+1)
          yfft(j)=yfft(i)
          yfft(j+1)=yfft(i+1)
          yfft(i)=tempr
          yfft(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*yfft(j)-sngl(wi)*yfft(j+1)
            tempi=sngl(wr)*yfft(j+1)+sngl(wi)*yfft(j)
            yfft(j)=yfft(i)-tempr
            yfft(j+1)=yfft(i+1)-tempi
            yfft(i)=yfft(i)+tempr
            yfft(i+1)=yfft(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

c-----------------------------------------------------------------------
      subroutine zkmeans_cluster(x0,y0,n,np,xscale,yscale,m,
     >xc,yc,vxc,vyc,bic)
c-----------------------------------------------------------------------
c
c      subroutine to find m clusters in a set of n pixels defined by x0 and y0. The x0 and y0 data are scaled by x/yscale. This is necessary because x0 and y0 may have very different ranges. For example if x is 0-0.0001 and y is 0-100 then the clustering will be dominated by the y axis. The algorithm is:
c      -set m trail cluster locations, equally spaced along the data
c      range diagonal
c      -asign each data point to it's nearest cluster
c      -calc the location of the cluster from the mean position of the
c      data points in that cluster
c      -re-asign the datapoints to the updated cluster locations
c      -continue until the asignment of points remains unchanged number
c       of itereations exceeds MAXITER
c    in:
c      x0(np)      real            x data
c      y0(np)      real            y data
c      n            int            number of pixels
c      np            int            array dimension
c      xscale      real            scale x data by this factor
c      yscale      real            scale y data by this factor
c      m            int            number of clusters
c    out:
c      xc/yc(np)      real            cluster positions (means)
c      vxc/vyc(np)      real            variance in x and y directions
c
c      Algorithm from:
c      Richards and Jia. Remote sensing digital image analysis, 1999,
c      springer (99RJ)
c
c-----------------------------------------------------------------------
c      N. Teanby      8-8-02      Original code
c-----------------------------------------------------------------------
      implicit none
      integer n,np,m,np_local,iter,maxiter
      real TINY
      parameter (np_local=50000,maxiter=100,TINY=1.e-32)
      real x0(np),y0(np),xc(np),yc(np),vxc(np),vyc(np),xscale,yscale,bic
      real x(np_local),y(np_local)
      real xmin,xmax,ymin,ymax,distmin,xsum,ysum,dist,term1,term2
      integer clusters(np_local),clustersOLD(np_local)
      integer i,j,nclust(np_local),nc
      logical clusters_stable

c  ** check np_local is big enough **
      if (n.gt.np_local) pause 'ERROR: zmigrating-means: n.gt.np_local'

c  ** scale the x and y values **
      do 11 i=1,n
         x(i)=x0(i)/xscale
         y(i)=y0(i)/yscale
11      continue

c  ** set initial cluster locations **
c  ** equally space along the diagonal from minimum x and y to the
c      maximum x and y (99RJ p227) **
c  ** find min and max of x and y **
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      do 1 i=1,n
         xmin=min(xmin,x(i))
         xmax=max(xmax,x(i))
         ymin=min(ymin,y(i))
         ymax=max(ymax,y(i))
1      continue
c  ** set the initial cluster positions **
      do 2 j=1,m
         xc(j) = xmin + real(j-1)*(xmax-xmin)/real(m-1)
         yc(j) = ymin + real(j-1)*(ymax-ymin)/real(m-1)
2      continue

      do 9 iter=1,maxiter
c      ** calculate the distance between the ith datapoint and the jth cluster **
c      ** assign each datapoint to the cluster which is closest **
         do 10 i=1,n
            distmin=sqrt((xmax-x(i))**2+(ymax-y(i))**2)
            do 20 j=1,m
               dist=sqrt((x(i)-xc(j))**2+(y(i)-yc(j))**2)
               if (dist.lt.distmin) then
                  distmin=dist
                  clusters(i)=j
               endif
20            continue
10         continue
c      ** if allocations are unchanged then exit the loop **
         clusters_stable=.true.
         do 25 i=1,n
           if (iter.gt.1) then
              if(clusters(i).ne.clustersOLD(i)) then
                 clusters_stable=.false.
              endif
           else
              clusters_stable=.false.
           endif
25         continue
         if (clusters_stable) then
            goto 99
         endif
c      ** recalculate the cluster centres **
         do 30 j=1,m
            nclust(j)=0
            xsum=0.
            ysum=0.
            do 40 i=1,n
               if (clusters(i).eq.j) then
                  xsum = xsum + x(i)
                  ysum = ysum + y(i)
                  nclust(j) = nclust(j) + 1
               endif
40            continue
            if (nclust(j).gt.0) then
               xc(j)=xsum/real(nclust(j))
               yc(j)=ysum/real(nclust(j))
            else
               xc(j)=0.
               yc(j)=0.
            endif
30         continue
c      ** copy cluster allocations **
         do 50 i=1,n
            clustersOLD(i)=clusters(i)
50         continue
9      continue
      print*,'ERROR: zmigrating-means:'
      pause 'maximum nuber of itererations reached, clusters not stable'
99      continue

c  ** remove the effect of scaling on the cluster centres **
      do 100 j=1,m
         xc(j)=xc(j)*xscale
         yc(j)=yc(j)*yscale
100      continue

c  ** calculate the variance of each cluster in the x and y directions **
      do 150 j=1,m
         nclust(j)=0
         xsum=0.
         ysum=0.
         do 300 i=1,n
               if (clusters(i).eq.j) then
                  xsum = xsum+ (x0(i)-xc(j))**2
                  ysum = ysum+ (y0(i)-yc(j))**2
               nclust(j)=nclust(j)+1
               endif
300         continue
         if (nclust(j).ne.0) then
            vxc(j)=xsum/real(nclust(j))
            vyc(j)=ysum/real(nclust(j))
         else
            vxc(j)=0.
            vyc(j)=0.
         endif
150      continue

c  ** calc Baysian Information Criterion **
      nc=0
      do 500 j=1,m
         if (nclust(j).ne.0) then
            nc=nc+1
         endif
500      continue
      bic=0.
      do 400 j=1,m
         if (nclust(j).ne.0) then
            term1=(vxc(j)/xscale + vyc(j)/yscale)/real(nclust(j))
            term2=nc*log(real(nclust(j)))/real(nclust(j))
            if (term1.gt.TINY) then
               bic = bic + log(term1) + term2
            endif
         endif
400      continue

      return
      end
c-----------------------------------------------------------------------
      subroutine zlag(x,y,n,np,lag,iwextra,xlag,ylag,noverlap)
c-----------------------------------------------------------------------
c
c      -LAG MUST BE +VE
c      -adds a lag to the x and y components
c      -positive lag => x has been shifted forward in time by lag points
c      -so y is the fast wave if lag is positive
c      -data from outside the analysis window are used such that the overlap
c      is constant (=analysis window length).
c
c      input data
c      <------analysis window------------>
c      |----------------------------------|-------|
c      1                             n-iwextra     n
c
c      output data
c      <------analysis window------------>
c      |----------------------------------|
c      1                                noverlap
c
c      variables
c      x(np)                  real      time series
c      y(np)                  real      time series
c      n                  int      number of points
c      np                  int      array dimension
c      lag                  int      lag (x shifted forward wrt y by lag points)
c      iwextra            int      number of extra points included in window
c
c      xlag(np)            real      lagged time series
c      ylag(np)            real      lagged time series
c      noverlap            int      length of lagged time series
c
c-----------------------------------------------------------------------
c      N. Teanby      30-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,lag,noverlap,noverlap_max,i,iwextra
      real x(np),y(np),xlag(np),ylag(np)

c  ** calc max number of overlapping points **
      noverlap_max = n - iwextra
      noverlap=0

c  ** lag the time series by lag, new time series are noverlap long
c      and the extra points have been omitted **
      if (lag.eq.0) then
         do i=1,noverlap_max
           xlag(i)=x(i)
           ylag(i)=y(i)
         enddo
         noverlap=noverlap_max
      else if (lag.gt.0) then
         do i=1,noverlap_max
            if (i+lag.le.n) then
               xlag(i)=x(i+lag)
               ylag(i)=y(i)
               noverlap=noverlap+1
            endif
         enddo
      else
         pause 'ERROR: zlag: negative lag not supported'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine zlinint(y,n,np,ninterp,yint)
c-----------------------------------------------------------------------
c
c      linear interpolate the series y from n points to ninterp points
c
c      variables
c      input:
c      y(np)                        real      series to interpolate
c      n                        int      number of points
c      np                        int      array dimension
c      ninterp                  int      number of points to interpolate to
c
c      yint(np)                  real      interpolated series
c
c-----------------------------------------------------------------------
c      N. Teanby      31-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer i,n,np,ninterp,index1,index2
      real y(np),yint(np)
      real x,dx

c  ** assume spacing of orignal series (y) = 1 **
c  ** spacing of interpolated series is then dx, where **
      dx = (n-1.)/(ninterp-1.)

c  ** interpolation **
      do 1 i=1,ninterp
c         print*,i
c      ** calc index of the two points used to do the interpolation **
         x      = 1. + ( real(i) - 1.) * dx
c         print*,x
         index1 = int( x )
         index2 = index1 + 1
c      ** interpolate **
         if (index1.eq.index2) then
            yint(i)=y(index1)
         else
            yint(i)=y(index1) + (y(index2)-y(index1))*(x-aint(x))
         endif
1      continue

      return
      end

c-----------------------------------------------------------------------
      subroutine zmultiply(a,n,np,scalar,axscalar)
c-----------------------------------------------------------------------
c
c      subroutine to multiply each element of 'a' by 'scalar' to
c      give 'axscalar'
c
c    in:
c      a(np)            real            array to be multiplied
c      n            int            number of data points
c      np            int            array dimension#
c      scalar      real            multiplication factor
c    out:
c      axscalar(np)real            scaled array
c
c-----------------------------------------------------------------------
c      N. Teanby      4-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i
      real a(np),scalar,axscalar(np)

      do 1 i=1,n
         axscalar(i) = a(i) * scalar
1      continue

      return
      end

c-----------------------------------------------------------------------
      subroutine zpackresults(wbeg,wend,tlag,dtlag,fast,dfast,n,npc,
     >dtlag_max,dfast_max)
c-----------------------------------------------------------------------
c
c      remove results which have errors over dtlag_max and dfast_max and
c      compress the arrays. results which only fail one of the criteria
c      are still accepted.
c
c    in/out:(modified on return)
c      wbeg(npc)      real            beginning of ith window
c      wend(npc)      real            end of ith window
c      fast(npc)      real            fast direction for ith window
c      dfast(npc)      real            s.d. of fast direction for ith window
c      tlag(npc)      real            tlag direction for ith window
c      dtlag(npc)      real            s.d. of tlag direction for ith window
c      n            int            number of data
c    in:
c      npc            int            array dimension
c      dtlag_max      real            max allowable error in lag time
c      dfast_max      real            max allowable error in fast direction
c
c-----------------------------------------------------------------------
c      n.teanby      20-8-02      original code
c-----------------------------------------------------------------------

      implicit none
      integer n,npc
      real wbeg(npc),wend(npc),fast(npc),dfast(npc)
      real tlag(npc),dtlag(npc)
      real dtlag_max,dfast_max
      integer i,j

c  ** remove bad results and compress arrays **
      j=0
      do i=1,n
         if ((dtlag(i).le.dtlag_max).or.(dfast(i).le.dfast_max)) then
            j=j+1
            wbeg(j)  = wbeg(i)
            wend(j)  = wend(i)
            tlag(j)  = tlag(i)
            dtlag(j) = dtlag(i)
            fast(j)  = fast(i)
            dfast(j) = dfast(i)
         endif
      enddo
      n=j

c  ** clean up end of array **
      do i=j+1,npc
         wbeg(i)  = 0.
         wend(i)  = 0.
         tlag(i)  = 0.
         dtlag(i) = 0.
         fast(i)  = 0.
         dfast(i) = 0.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zreaddataSAC(file1,file2,lu,force_spick,
     >x,y,n,np,b,delta,as,fs,ppick,spick,phi,theta,
     >ev_x,ev_y,ev_z,ev_time)
c-----------------------------------------------------------------------
c      read in two binary sac files containing two components of a seismogram
c
c      force_spick      logical      .true. ERROR if spick not defined
c                        (use in zass.f, where window is defined using spick)
c                              .false. no error if spick not defined
c                        (use in split.f, where window is defined)
c
c      x/y(np)      real      data in files 1 and 2
c      n            int      number of data points
c      np            int      array dimensions
c      b            real      start of time series
c      delta            real      sampling interval
c      as            real      beginning of s-wave window (hand pick)
c      fs            real      end of s-wave window (hand pick)
c      ppick            real      p-wave onset (seconds) read from T4 header variable
c      spick            real      s-wave onset (seconds) read from T5 header variable
c      phi            real      rotation angle of frame (deg clock from N)
c      theta            real      rotation angle of frame (angle between z and z')
c      ev_x/y/z      real      x/y/z position of event read from sac header in metres
c                        (stored in RESP4/5/6, not ideal but no x y z loc in
c                        SAC header)
c      ev_time      real      time of event in days
c
c-----------------------------------------------------------------------
c      N. Teanby      20-8-02      original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,lu
      real x(np),y(np),b,delta,as,fs,ppick,spick,phi,theta
      real ev_x,ev_y,ev_z,ev_time
      logical force_spick
c  ** header for sac file
      real h1(14,5)
      integer h2(8,5)
      character h3a*8, h3b*16
      character*8 h4(7,3)
c  **
      character*50 file1,file2

c  ** read in data (SAC binary file) **
      call zbreadsac(file1,lu,np,h1,h2,h3a,h3b,h4,x,n)
      call zbreadsac(file2,lu,np,h1,h2,h3a,h3b,h4,y,n)

c  ** extract header info **
      delta  = h1(1,1)
      b      = h1(2,1)
      as     = h1(3,3)
      fs     = h1(3,4)
      ppick  = h1(3,5)
      spick  = h1(4,1)
      phi       = h1(5,2)
      theta       = h1(5,4)
      ev_x   = h1(6,1)
      ev_y   = h1(6,2)
      ev_z   = h1(6,3)
      ev_time= real(h2(1,2)) + real(h2(1,3))/24. + real(h2(1,4))/1440. +
     > real(h2(1,5))/86400. + real(h2(2,1))/86400000.

c  ** make sure the essential quantities are defined **
      if (delta.lt.0.) pause 'ERROR: zreaddataSAC: delta= -ve'
      if (b    .lt.0.) pause 'ERROR: zreaddataSAC: b    = -ve'
      if (force_spick) then
        if (spick.lt.0.) pause 'ERROR: zreaddataSAC: spick= -ve'
      endif

c  ** check other non-essential quantities are defined **
      if (nint(phi).eq.-12345) then
         phi=0.
         print*,'WARNING: zreaddataSAC: phi undefined - set to 0'
      endif
      if (nint(theta).eq.-12345) then
         theta=0.
         print*,'WARNING: zreaddataSAC: theta undefined - set to 0'
      endif
      if (nint(ppick).eq.-12345) then
         print*,'WARNING: zreaddataSAC: ppick = undefined'
      endif
      if (nint(as).eq.-12345) then
         print*,'WARNING: zreaddataSAC: as = undefined'
      endif
      if (nint(fs).eq.-12345) then
         print*,'WARNING: zreaddataSAC: fs = undefined'
      endif
      if (h2(2,1).eq.-12345) then
         ev_time =0.
         print*,'WARNING: zreaddataSAC: ev_time undefined, set to 0'
      endif
      if ((nint(ev_x).eq.-12345).or.(nint(ev_y).eq.-12345).or.
     >      (nint(ev_z).eq.-12345)) then
         ev_x=0.
         ev_y=0.
         ev_z=0.
         print*,'WARNING: zreaddataSAC: ev_x/y/z undefined, set to 0'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine zreadsac(file,lu,np,h1,h2,h3a,h3b,h4,d1,npts)
c-----------------------------------------------------------------------
c
c      read in a sac ascii data file
c        read in header info (in 4 parts)
c        read in data part1
c        (may need to add a data part2 if unevenly spaced data is ever used)
c
c      variables:
c input: file      char*50            file to read in
c         lu            int                  logical unit to read file to
c         np            int                  dimension of d1 array
c output:h1            real(14,5)            header part1
c         h2            int(8,5)            header part2
c         h3a      char*8            header part3a
c         h3b      char*16            header part3b
c         h4            char*8 (7,3)      header part4
c         d1            real(np)            data points
c         npts      int                  number of data points (read from h2(2,5))
c
c-----------------------------------------------------------------------
c  modifications:
c      11-07-01      N. Teanby      Original code
c-----------------------------------------------------------------------

      implicit none
      integer lu,np
      real h1(14,5)
      integer h2(8,5)
      character h3a*8, h3b*16
      character*8 h4(7,3)
      real d1(np)
      integer i,npts
      character*50 file

c  ** open data file **
      open(lu,file=file,status='old')

c  ** read header info **
      do 1 i=1,14
         read(lu,1001) h1(i,1),h1(i,2),h1(i,3),h1(i,4),h1(i,5)
ccc         write(6,1001) h1(i,1),h1(i,2),h1(i,3),h1(i,4),h1(i,5)
1      continue
      do 2 i=1,8
         read(lu,1002) h2(i,1),h2(i,2),h2(i,3),h2(i,4),h2(i,5)
ccc         write(6,1002) h2(i,1),h2(i,2),h2(i,3),h2(i,4),h2(i,5)
2      continue
      read(lu,1003) h3a,h3b
ccc      write(6,1003) h3a,h3b
      do 4 i=1,7
         read(lu,1004) h4(i,1),h4(i,2),h4(i,3)
ccc         write(6,1004) h4(i,1),h4(i,2),h4(i,3)
4      continue

c  ** read npts from header **
      npts=h2(2,5)
c  ** check npts does not exceed array dimension **
      if (npts.gt.np) then
         print*,'WARNING:zsacread'
         print*,'WARNING:  npts.gt.np, data missed out'
      endif

c  ** read in data **
      do 100 i=1,5*int(npts/5),5
         read(lu,1010) d1(i),d1(i+1),d1(i+2),d1(i+3),d1(i+4)
ccc         write(6,1010),d1(i),d1(i+1),d1(i+2),d1(i+3),d1(i+4)
100      continue
c  ** check that npts is a multiple of 5 **
      if (int(npts/5).ne.npts/5) then
         print*,'WARNING:zsacread'
         print*,'WARNING:  npts not multiple of 5'
         print*,'WARNING:  last few (<5) points may be missed off'
      endif

      close(lu)
      return

c **      format statements **
1001      format(5g15.7)
1002      format(5i10)
1003      format(a8,a16)
1004      format(3a8)
1010      format(5g15.7)

      end

c-----------------------------------------------------------------------
      subroutine zresample(a,n,np,f,ar,nr)
c-----------------------------------------------------------------------
c
c      subroutine to reduce the sampling of 'a' by a factor of f
c
c    in:
c      a(np)            real            array to be resampled
c      n            int            number of data points
c      np            int            array dimension
c      f            int            resampling factor
c    out:
c      ar(np)      real            resampled array
c      nr            int            number of resampled points
c
c-----------------------------------------------------------------------
c      N. Teanby      6-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i,nr,f
      real a(np),ar(np)

      if (f.lt.1) then
         pause 'ERROR: zresample: resampling factor f.lt.1'
      endif

c  ** number of resampled points **
      nr = int(real(n-1)/real(f))+1

      do 1 i=1,nr
         ar(i) = a(1+(i-1)*f)
1      continue

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine zrotabc2enz(phi,theta,angle_abc,dangle_abc,
     >angle_strike_enz,dangle_strike_enz,angle_dip_enz,dangle_dip_enz)
c-----------------------------------------------------------------------
c
c      -find strike and dip of a plane
c
c      -the original e n z coords have been rotated into a b c (in frame of ray)
c      -the angles az and polaz are measured in the a b c frame
c      -a b c is rotated by phi and theta relative to e n z
c      -this prog asumes that angle abc defines a plane parralell to c axis
c      -the normal to this plane rotated into the e n z frame and used
c       to define the plane in the enz frame (plane goes through origin in each
c       case)
c      -angle enz is the angle, clockwise from north, of the intersection of the
c       plane with the en horizontal plane
c
c-----------------------------------------------------------------------
c      n. teanby      31-08-01      original code
c-----------------------------------------------------------------------

      implicit none
      real deg2rad,fdstr,fddip
      parameter (deg2rad=0.0174533)
      real phi,theta,angle_abc,dangle_abc
      real angle_strike_enz,dangle_strike_enz
      real angle_dip_enz,dangle_dip_enz
      real phi_rad,theta_rad,angle_abc_rad,nx,ny,nz

c  ** convert angles to radians **
      phi_rad       = phi       * deg2rad
      theta_rad     = theta     * deg2rad
      angle_abc_rad = angle_abc * deg2rad

c  ** calculate the normal to the plane in enz coords **
      nx =   cos(phi_rad)*cos(angle_abc_rad)
     >     - sin(phi_rad)*cos(theta_rad)*sin(angle_abc_rad)
      ny = - sin(phi_rad)*cos(angle_abc_rad)
     >     - cos(phi_rad)*cos(theta_rad)*sin(angle_abc_rad)
      nz =   sin(theta_rad)*sin(angle_abc_rad)

c  ** calc strike and dip **
      call zstrikedip(nx,ny,nz,angle_strike_enz,angle_dip_enz)

c  ** get strike in range 0-180 degrees (initially -180 to 180)**
      if (angle_strike_enz.lt.0.) then
         angle_strike_enz = angle_strike_enz + 180.
      endif

c  ** error on strike and dip of angle in enz **
      dangle_strike_enz =
     >( fdstr(theta,angle_abc, dangle_abc) +
     >  fdstr(theta,angle_abc,-dangle_abc) )/2.
      dangle_dip_enz   =
     >( fddip(theta,angle_abc, dangle_abc) +
     >  fddip(theta,angle_abc,-dangle_abc) )/2.

      return
      end

c-----------------------------------------------------------------------
      subroutine zstrikedip(nx,ny,nz,strike,dip)
c-----------------------------------------------------------------------
c
c      -find strike and dip of a plane defined by the plane normal
c      (nx,ny,nz)
c
c      output strike and dip in degrees
c-----------------------------------------------------------------------
c      n. teanby      18-09-01      original code
c-----------------------------------------------------------------------

      implicit none
      real pi
      parameter (pi=3.141592654)
      real nx,ny,nz,strike,dip

c  ** find dip (rem vector has magnitude = 1, use cos = base/hypot **
      dip = asin( sqrt(nx**2 + ny**2))

c  ** find strike **
      if (nz.ge.0.) then
         strike = atan2(-ny,nx)
      else
         strike = atan2(ny,-nx)
      endif

c  ** convert to degrees **
      dip    = dip    * 180./pi
      strike = strike * 180./pi

      return
      end

c-----------------------------------------------------------------------
      function fdstr(theta,az,daz)
c-----------------------------------------------------------------------
c
c      consider a plane in abc frame defined by angle az
c      consider a second plane in abc frame defined by az+daz
c
c      this function find the angle between the strikes of the two
c      planes in the enz frame.
c
c      expression calculated by defining the plane normals in abc
c      rotating these normals to enz and finding the angle between
c      the horizontal components using the dot product
c
c in:      theta            real      rotation angle theta in degrees
c      az            real      azimuth of plane 1 in abc frame in degrees
c      daz            real      difference between azimuths of plane 1 and 2
c                        in abc frame in degrees
c      RETURN      real      difference between strike of two planes in
c                        enz frame in degrees
c
c-----------------------------------------------------------------------
c      n. teanby      12-12-01      original code
c-----------------------------------------------------------------------
      implicit none

      real pi
      parameter (pi=3.141592654)
      real theta,az,daz,fdstr
      real theta_rad, az_rad, daz_rad
      real dotprod,mod

      theta_rad      = theta      *pi/180.
      az_rad      = az            *pi/180.
      daz_rad      = daz            *pi/180.

      dotprod =       cos(az_rad)*cos(az_rad+daz_rad) +
     >            cos(theta_rad)**2 * sin(az_rad)*sin(az_rad+daz_rad)

      mod = sqrt((cos(az_rad)**2 + cos(theta_rad)**2 * sin(az_rad)**2)*
     >(cos(az_rad+daz_rad)**2+cos(theta_rad)**2*sin(az_rad+daz_rad)**2))

      fdstr    =      acos(dotprod/mod)*180./pi

      return
      end

c-----------------------------------------------------------------------
      function fddip(theta,az,daz)
c-----------------------------------------------------------------------
c
c      consider a plane in abc frame defined by angle az
c      consider a second plane in abc frame defined by az+daz
c
c      this function find the angle between the dips of the two
c      planes in the enz frame.
c
c      expression calculated by defining the vertical comp of the
c      plane normals in abcthis  to enz and finding the angle between
c      the vertical components
c
c in:      theta            real      rotation angle theta in degrees
c      az            real      azimuth of plane 1 in abc frame in degrees
c      daz            real      difference between azimuths of plane 1 and 2
c                        in abc frame in degrees
c      RETURN      real      difference between dips of two planes in
c                        enz frame in degrees
c
c-----------------------------------------------------------------------
c      n. teanby      12-12-01      original code
c-----------------------------------------------------------------------
      implicit none

      real pi
      parameter (pi=3.141592654)
      real theta,az,daz,fddip
      real theta_rad, az_rad, daz_rad
      real nz,nz0

      theta_rad      = theta      *pi/180.
      az_rad      = az            *pi/180.
      daz_rad      = daz            *pi/180.

c  ** calc vertical comp of plane normal in enz frame **
      nz0 =   sin(theta_rad)*sin(az_rad)
      nz  =   sin(theta_rad)*sin(az_rad + az_rad)

      fddip    =      abs(acos(nz0)-acos(nz))*180./pi

      return
      end
c-----------------------------------------------------------------------
      subroutine zrotate2d(x,y,n,np,angle,xrot,yrot)
c-----------------------------------------------------------------------
c
c      rotate a 2 seismogram (x and y comps) by angle.
c      angle is in degrees measured clockwise
c
c      variables
c      x(np)                        real      time series
c      y(np)                        real      time series
c      n                        int      number of points
c      np                        int      array dimension
c      angle                        real      rotation angle (clockwise from North (y))
c
c      xrot(np)                  real      rotated time series
c      yrot(np)                  real      rotated time series
c
c-----------------------------------------------------------------------
c      N. Teanby      30-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i
      real pi
      parameter (pi=3.141592654)
      real x(np),y(np),xrot(np),yrot(np)
      real angle,angler

c  ** convert to radians **
      angler=angle*pi/180.

c  ** rotate time series x and y by angle **
      do 1 i=1,n
         xrot(i) = x(i)*cos(angler) - y(i)*sin(angler)
         yrot(i) = x(i)*sin(angler) + y(i)*cos(angler)
1      continue

      return
      end

c-----------------------------------------------------------------------
      subroutine zselect_cluster(dx0,dy0,vxc0,vyc0,n,
     >xscale,yscale,cluster,nmin,k,
     >kbest)
c-----------------------------------------------------------------------
c
c      subroutine to find the best cluster.
c
c      best cluster has:      nmin or greater points; and the lowest overall variance
c
c    in:
c      dx0/dy0(npc)      real      standard deviation of tlag and fast measurements
c      vxc0/vyc0(npc)      real      within cluster variance
c      n                  int      number of data points
c      x/yscale            real      scale/standardisation factors for x0/y0 data
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c      nmin                  int      min no. data points for an acceptable cluster
c      k                  int      number of clusters
c    out:
c      kbest                  int      best cluster
c    other:
c      dx/dy(npc)            real      scaled data standard deviation
c      vxc/vyc(npc,npc)      real      scaled cluster variance of ith cluster when k=j
c      nc(npc)            int      number of data points in ith cluster
c      var_data(npc)      real      average variance of data within ith cluster
c      var_cluster(npc)      real      within cluster varaince (=vxc+vyc)
c      var_overall(npc)      real      overall variance = max(var_cluster,var_data)
c-----------------------------------------------------------------------
c      n.teanby      15-8-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,j,k,kbest
      real dx0(npc),dy0(npc),vxc0(npc),vyc0(npc)
      real dx(npc),dy(npc),vxc(npc,npc),vyc(npc,npc)
      real xscale,yscale
      integer nc(npc),cluster(npc,npc)
      real var_overall(npc),var_data(npc),var_cluster(npc),var_min
      integer nmin
      logical first_pass

      kbest=0

c  ** scale the errors on the data **
      call zmultiply(dx0,n,npc,1./xscale,dx)
      call zmultiply(dy0,n,npc,1./yscale,dy)

c  ** calc composite variance of data in clusters **
      call zcluster_variance(dx,dy,n,npc,cluster,k,var_data)

c  ** calc variance of points in cluster about the cluster mean **
      do j=1,k
         var_cluster(j) = vxc(j,k)**2 + vyc(j,k)**2
      enddo

c  ** calc number of data, nc, in each of the k clusters **
      call zcluster_number(cluster,n,npc,k,nc)

c  ** set variance to be maximum of var_data, and var_cluster **
      do j=1,k
         var_overall(j) = max(var_cluster(j),var_data(j))
         print*,'cluster no.=',j,' error = ',var_overall(j),' n=',nc(j)
      enddo

c  ** find best cluster (lowest error and n>=nmin) **
      first_pass=.true.
      do j=1,k
         if(nc(j).ge.nmin) then
            print*,'cluster ',j,' passes n>',nmin
            if (first_pass) then
               first_pass = .false.
               var_min    = var_overall(j)
               kbest      = j
            else if(var_overall(j).lt.var_min) then
               var_min = var_overall(j)
               kbest   = j
            endif
         endif
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine zcluster_variance(dx,dy,n,npc,cluster,k,var_data)
c-----------------------------------------------------------------------
c
c      variance =        1
c                 --------------
c                 SUM( 1/var_i )
c
c      using this definition of the variance (instead of the average
c      variance) means that the overall variance is dominated by the
c      smaller individual variances.
c
c-----------------------------------------------------------------------
c      n.teanby      15-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,npc
      real dx(npc),dy(npc),var_data(npc),TINY
      parameter (TINY=1.e-20)
      integer k,cluster(npc,npc)
      real xsum,ysum
      integer nc,i,j

c  ** find variance of data in jth cluster **
      do j=1,k
         xsum=0.
         ysum=0.
         nc=0
         do i=1,n
            if (cluster(i,k).eq.j) then
c            ** don't include zero variance points **
               if ((dx(i).gt.TINY).and.(dy(i).gt.TINY)) then
                  xsum = xsum + 1./dx(i)**2
                  ysum = ysum + 1./dy(i)**2
                  nc   = nc+1
               endif
            endif
         enddo
         var_data(j) = 1./xsum + 1./ysum
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine zcluster_number(cluster,n,npc,k,nc)
c-----------------------------------------------------------------------
c
c      count the number of data points in each cluster
c
c-----------------------------------------------------------------------
c      n.teanby      15-8-02      original code
c-----------------------------------------------------------------------
      implicit none
      integer n,npc,k
      integer cluster(npc,npc),nc(npc)
      integer i,j

c  ** calc number of data in jth cluster **
      do j=1,k
         nc(j)=0
         do i=1,n
            if (cluster(i,k).eq.j) then
               nc(j) = nc(j) + 1
            endif
         enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine zselect_measurement(dx0,dy0,n,
     >wbeg,wend,delta,spick,xscale,yscale,cluster,k,kbest,
     >ibest)
c-----------------------------------------------------------------------
c
c      subroutine to find the best measurement from
c      within the best cluster.
c
c    in:
c      dx0/dy0(npc)      real      standard deviation of tlag and fast measurements
c      n                  int      number of data points
c      delta            real            sampling interval (s) read from sac header
c      spick            real            s-wave pick read from sac header
c      x/yscale            real      scale/standardisation factors for x0/y0 data
c      cluster(npc,npc)      int      assignment of datapoints. eg cluster(3,17) is the
c                              number of the cluster that the 3rd datapoint
c                              is in for 17 clusters
c      k                  int      optimum number of clusters = max(k1,k2)
c      kbest                  int      best cluster
c    out:
c      ibest                  int      index of best measurement within best cluster
c    other:
c      dx/dy(npc)            real      scaled data standard deviation
c-----------------------------------------------------------------------
c      n.teanby      15-8-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,k,kbest,ibest
      real dx0(npc),dy0(npc),dx(npc),dy(npc)
      real xscale,yscale,wbeg(npc),wend(npc),delta,spick
      integer cluster(npc,npc)
      logical first_pass
      integer i
      real err,err_min

      ibest=0

      if (kbest.ne.0) then
c        ** scale the data **
         call zmultiply(dx0,n,npc,1./xscale,dx)
         call zmultiply(dy0,n,npc,1./yscale,dy)
c        ** find best splitting measurement from the best cluster **
         first_pass=.true.
         do i=1,n
            if (cluster(i,k).eq.kbest) then
               err=sqrt(dy(i)**2 + dx(i)**2)/((wend(i)-spick)/delta)
c               err=sqrt(dy(i)**2 + dx(i)**2)
               if (first_pass) then
                  first_pass = .false.
                  ibest      = i
                  err_min    = err
               else if (err.le.err_min) then
                  ibest   = i
                  err_min = err
               endif
            endif
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine zsourcepol(fast,lambda1,lambda2,vec1,vec2,spol,dspol)
c-----------------------------------------------------------------------
c
c      -calculate source polarisation and error
c      -error is the MAD (max angular deviation) and is the angle
c      subtended by the min and max eigenvalues (see kirschvink 1980, PCA)
c
c      variables
c      fast                  real      fast direction (deg clock from N)
c      lambda1/2            real      largest/smallest eignevalue
c      vec1/2            real      largest/smallest eignevector
c      spol                  real      polarisation direction (deg clock from N)
c      dspol                  real      error (MAD angle)
c
c-----------------------------------------------------------------------
c      N. Teanby      31-7-02      Original code
c      J. Wookey      6-10-08      Modified for sensible angle
c-----------------------------------------------------------------------

      implicit none
      real fast,lambda1,lambda2,vec1(2),vec2(2),spol,dspol,pi
      parameter (pi=3.141592654)

c  ** polarization in azimuth degrees(clockwise from north)
c      note, eigenvectors are in rotated coord frame and must add
c      fast to the angle. **
      spol = fast + atan2(vec1(1),vec1(2))*180./pi

c  ** calc MAD angle as a measure of the error **
      dspol = atan(lambda2/lambda1)*180./pi

c  ** recast into 0-360
      do ! forever
         if (spol.ge.0 .and. spol.lt.360) exit  
         if (spol.ge.360.) spol = spol - 360.0
         if (spol.lt.0) spol = spol + 360.0   
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine zsplint(y,n,f,yinterp,ninterp)
c-----------------------------------------------------------------------
c      interpolate y by integer factor of f using NR cubic B-spline routines
c
c      variables
c    in:
c      y(np)            real      data series to interpolate
c      np            int      array dimension (read from SIZE_np.h at compile time)
c      n            int      number of data points
c      f            int      interpolation factor
c    out:
c      yinterp(np)      real      interpolated y
c      ninterp      int      number of interpolated points = f(n-1)+1
c    local:
c      y2(np)      int      second derivatives of splines
c
c-----------------------------------------------------------------------
c      N. Teanby      21-8-02      Original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,f,i,ninterp
      real y(np),yinterp(np),y2(np),x

c  ** calc number of interpolated datapoints **
      ninterp=f*(n-1)+1
      if (ninterp.gt.np) then
         print*,'ninterp=',ninterp
         pause 'ERROR: zsplint: resampling creates too many datapoints'
      endif

c  ** find second derivatives of interpolating splines **
      call spline(y,n,y2)

c  ** do interpolation **
      do 1,i=1,ninterp
         x=1+real(i-1)/real(f)
         call splint(y,y2,n,np,x,yinterp(i))
1      continue

      return
      end

      SUBROUTINE spline(y,n,y2)
c-----------------------------------------------------------------------
c  ** modified for natural spline only and unit x spacing **
c   y = data to interpolate
c   n = number of data points
c   np = array dimensions
c   y2 = 2nd derivatives of splines
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      INTEGER n
      REAL y(np),y2(np)
      INTEGER i,k
      REAL p,qn,sig,un,u(np)

      y2(1)=0.
      u(1)=0.

      do 11 i=2,n-1
        sig=1./2.
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))-(y(i)-y(i-1)))/2.-sig*u(i-1))/p
11    continue
      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

c-----------------------------------------------------------------------
c  ** modified for unit x spacing **
c   ya = data to interpolate
c   n = number of data points
c   np = array dimensions
c   y2a = 2nd derivatives of splines (from spline)
c   x = x value to interpolate at
c   y = interpolated value
      SUBROUTINE splint(ya,y2a,n,np,x,y)
      implicit none
      INTEGER n,np
      REAL x,y,y2a(np),ya(np)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(real(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=real(khi-klo)
      if (h.eq.0.) pause 'bad x input in splint'
      a=(khi-x)/h
      b=(x-real(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

c-----------------------------------------------------------------------
      subroutine zsplit_ini(inifile,lu,
     >      ext1,ext2,fast_scale,tlag_scale,OPT_outfiles)
c-----------------------------------------------------------------------
c      read in fast_scale and tlag_scale from the split.ini file
c
c      in:
c       inifile      char*50      initialisation file
c      out:
c       ext1/2      char*50      extention of the two components to do analysis on
c       fast_scale      real            grid search scale for fast direction
c       tlag_scale      real            grid search scale for tlag
c       OPT_outfiles logical      true if write outfiles for gmt plots
c
c-----------------------------------------------------------------------
c      N. Teanby      12-5-03      Original code
c-----------------------------------------------------------------------

      implicit none
      integer lu
      real fast_scale,tlag_scale
      logical OPT_outfiles
      character*50 inifile,ext1,ext2

      open(lu,file=inifile,status='old')
      read(lu,*) ext1
      read(lu,*) ext2
      read(lu,*) fast_scale
      read(lu,*) tlag_scale
      read(lu,*) OPT_outfiles
      close(lu)

      return
      end
c-----------------------------------------------------------------------
      subroutine zwindow(x,y,n,np,wbeg,wend,xwindow,ywindow,nwindow)
c-----------------------------------------------------------------------
c
c      window time series x and y between index wbeg and wend
c
c      variables
c      x(np)                        real      time series
c      y(np)                        real      time series
c      n                        int      number of points
c      np                        int      array dimension
c      wbeg                        int      begin window
c      wend                        int      end window
c
c      xwindow(np)                  real      windowed time series
c      ywindow(np)                  real      windowed time series
c      nwindow                  int      length of windowed times series
c
c-----------------------------------------------------------------------
c      N. Teanby      30-7-02      Original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i,nwindow,wbeg,wend
      real x(np),y(np),xwindow(np),ywindow(np)

      nwindow = wend - wbeg + 1

c  ** window time series **
      do 1 i=1,nwindow
         xwindow(i) = x(i+wbeg-1)
         ywindow(i) = y(i+wbeg-1)
1      continue

      return
      end

c-----------------------------------------------------------------------
      subroutine zwrite1(file,lu,fmt,x,n,np)
c-----------------------------------------------------------------------
c      write out a 1 column datafile
c-----------------------------------------------------------------------
c      n. teanby      22/10/01      original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i,lu
      real x(np)
      character*50 file
      character*50 fmt

      open(lu,file=file,status='unknown')
      do i = 1,n
         write(lu,fmt) x(i)
      enddo
      close(lu)

      return
      end
c-----------------------------------------------------------------------
      subroutine zwrite2(file,lu,fmt,x,y,n,np)
c-----------------------------------------------------------------------
c      write out a 2 column datafile
c-----------------------------------------------------------------------
c      n. teanby      22/10/01      original code
c-----------------------------------------------------------------------

      implicit none
      integer n,np,i,lu
      real x(np),y(np)
      character*50 file
      character*50 fmt

      open(lu,file=file,status='unknown')
      do i = 1,n
         write(lu,fmt) x(i),y(i)
      enddo
      close(lu)

      return
      end
c-----------------------------------------------------------------------
      subroutine zwrite_logfile(event,lu,file_log,
     >      ev_x,ev_y,ev_z,ev_time,
     >      wbeg_best,wend_best,tlag_best,dtlag_best,
     >      fast_best,dfast_best,spol_best,dspol_best,
     >      fastSTR,dfastSTR,fastDIP,dfastDIP,
     >      spolSTR,dspolSTR,spolDIP,dspolDIP)
c-----------------------------------------------------------------------
c      write out log file
c-----------------------------------------------------------------------
c      n.teanby      9-9-02      original code
c-----------------------------------------------------------------------

      implicit none
      integer lu
      real wbeg_best,wend_best
      real tlag_best,dtlag_best
      real fast_best,dfast_best,spol_best,dspol_best
      real fastSTR,dfastSTR,fastDIP,dfastDIP
      real spolSTR,dspolSTR,spolDIP,dspolDIP
      real ev_x,ev_y,ev_z,ev_time
      character*50 file_log,event

c  ** write log file **
      open(lu,file=file_log,status='unknown')
      write(lu,1021)
      write(lu,1020) event,ev_x,ev_y,ev_z,ev_time,
     >   wbeg_best,wend_best,tlag_best,dtlag_best,
     >   fast_best,dfast_best,spol_best,dspol_best,
     >   fastSTR,dfastSTR,fastDIP,dfastDIP,
     >   spolSTR,dspolSTR,spolDIP,dspolDIP
      close(lu)

c  ** formats **
1020      format(a50,2x,e15.7,2x,e15.7,2x,e15.7,2x,e15.7,2x,
     >   f10.6,2x,f10.6,2x,f10.6,2x,f10.6,2x,
     >   f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,
     >   f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,
     >   f8.3,2x,f8.3,2x,f8.3,2x,f8.3)
1021      format('event',61x,'x',16x,'y',16x,'z',13x,'time',3x,
     >   'wbeg_best',3x,'wend_best',8x,
     >   'tlag',7x,'dtlag',6x,
     >   'fast ',4x,'dfast',7x,'spol',5x,'dspol',3x,
     >   'fastSTR',2x,'dfastSTR',3x,'fastDIP',2x,'dfastDIP',3x,
     >   'spolSTR',2x,'dspolSTR',3x,'spolDIP',2x,'dspolDIP')

      return
      end
c-----------------------------------------------------------------------
      subroutine zwrite_outfiles(event,lu,tlag_scale,fast_scale,
     >x0,y0,n,b,delta,ppick,spick,
     >wbeg,wend,
     >fast,tlag,spol,error_int)
c-----------------------------------------------------------------------
c      write output files
c-----------------------------------------------------------------------
c      N. Teanby      6-8-02      original code
c-----------------------------------------------------------------------
      use array_sizes ! use the array_sizes modules
C-----------------------------------------------------------------------
      implicit none
      integer n,ninterp,noverlap,noverlap_resamp,nw,lu
      real wbeg,wend,delta,b,fast,tlag,spol,ppick,spick
      real error_int(np1,np2int),tlag_scale,fast_scale
c  **
      character*50 event
      integer f,iwbeg,iwend,ilag,itlag_step
      real x(np),y(np),x0(np),y0(np)
      real xinterp(np),yinterp(np),xw(np),yw(np)
      real xrot(np),yrot(np),xlag(np),ylag(np)
      real time(np),timew(np),sto(np),sro(np),stc(np),src(np)
      real stc_resamp(np),src_resamp(np)
      real ssw(np),sfw(np),sscw(np),sfcw(np)
      real sswn(np),sfwn(np),sscwn(np),sfcwn(np)
      real xwn(np),ywn(np),xcw(np),ycw(np),xcwn(np),ycwn(np)
c  ** filenames **
      character*50 file_sto,file_sro,file_stc,file_src,file_ss,file_sf
      character*50 file_ssc,file_sfc,file_pm,file_pmc,file_error
      character*50 file_gmt
c  ** other **
      real test_tlag,test_fast,norm_pm,cov(2,2),ftick
      real max_xw,max_yw,max_xcw,max_ycw
      real max_ssw,max_sfw,max_sscw,max_sfcw
      real max_sto,max_sro,max_stc,max_src,max_rt
      character*50 fmt2
      character*50 fstrcat,ext
      external fstrcat
      integer i,j
c  ** added by James Wookey + Jen Caddick ***
      real gmt_trace_start,gmt_trace_end

c  ** set filenames **
      ext='.sto'
      file_sto   = fstrcat(event,ext)
      ext='.sro'
      file_sro   = fstrcat(event,ext)
      ext='.stc'
      file_stc   = fstrcat(event,ext)
      ext='.src'
      file_src   = fstrcat(event,ext)
      ext='.ss'
      file_ss    = fstrcat(event,ext)
      ext='.sf'
      file_sf    = fstrcat(event,ext)
      ext='.ssc'
      file_ssc   = fstrcat(event,ext)
      ext='.sfc'
      file_sfc   = fstrcat(event,ext)
      ext='.pm'
      file_pm    = fstrcat(event,ext)
      ext='.pmc'
      file_pmc   = fstrcat(event,ext)
      ext='.error'
      file_error = fstrcat(event,ext)
      ext='.gmt'
      file_gmt   = fstrcat(event,ext)

c  ** format for 2 column output **
      fmt2='(e15.7,2x,e15.7)'

c--------------------- DATA ---------------------
c  ** calculateinterpolation factor for error surface from np2 and np2int **
c  ** this factor is used to interpolate the error surface in the tlag dirn **
      f = (np2int-1)/(np2-1)
c  ** f needs to be an integer for zsplint to work **
      if (mod(np2int-1,np2-1).ne.0) then
         pause 'ERROR: zwrite_outfiles: f not a whole number'
      endif
c  ** calc itlag_step from tlag_scale **
c  ** itlag_step is the grid spacing in tlag for the grid search **
c  ** it must be an integer greater than 1 **
      itlag_step = nint( tlag_scale/(real(np2-1)*delta) )

c  ** detrend the data **
      call zdetrend(x0,n,np,x)
      call zdetrend(y0,n,np,y)
c  ** set up an array with the times in **
      do i=1,n
         time(i) = b + real(i-1)*delta
      enddo
c  ** interpolate the time series by the factor f **
      call zsplint(x,n,f,xinterp,ninterp)
      call zsplint(y,n,f,yinterp,ninterp)

c------------------- OUTFILES -------------------
c*****UNCORRECTED RADIAL AND TRANSVERSE **
c  ** rotate into the radial and transverse directions **
c  ** RADIAL     [sro] = parallel to the source polarisation (spol) **
c  ** TRANSVERSE [sto] = perpendicular to the source polarisation (spol) **
      call zrotate2d(x,y,n,np,spol,sto,sro)
      call zwrite2(file_sto,lu,fmt2,time,sto,n,np)
      call zwrite2(file_sro,lu,fmt2,time,sro,n,np)

c*****CORRECTED RADIAL AND TRANSVERSE **
c  ** rotate into fast and slow direction and apply lag, then rotate
c      into radial and transverse directions **
c  ** because tlag may not be an integer number of samples, interpolation
c      is required **
c  ** calc the lag in terms of the interpolated index **
      ilag = 1 + nint(tlag*real(f)/delta)
c  ** rotate and lag (positive lag means that both time series still begin
c      at b, however the time series will be tlag shorter than the
c      uncorrected ones**
      call zrotate2d(xinterp,yinterp,ninterp,np,fast,xrot,yrot)
      call zlag(xrot,yrot,ninterp,np,ilag,0,xlag,ylag,noverlap)
c  ** rotate into spol dirn (radial and transverse) **
      call zrotate2d(xlag,ylag,noverlap,np,spol-fast,stc,src)
c  ** note these time series are interpolated by the factor f **
      call zresample(stc,noverlap,np,f,stc_resamp,noverlap_resamp)
      call zresample(src,noverlap,np,f,src_resamp,noverlap_resamp)
      call zwrite2(file_stc,lu,fmt2,time,stc_resamp,noverlap_resamp,np)
      call zwrite2(file_src,lu,fmt2,time,src_resamp,noverlap_resamp,np)
c  ** find maximum of r/t comps for the perposes of plotting the graphs **
      call zabs_max(sto,noverlap,np,max_sto)
      call zabs_max(sro,noverlap,np,max_sro)
      call zabs_max(stc,noverlap,np,max_stc)
      call zabs_max(src,noverlap,np,max_src)
      max_rt=max(max_sto,max_sro,max_stc,max_src)

c*****FAST AND SLOW COMPONENTS (WINDOWED)**
c  ** window the interpolated data and set up windowed time matrix **
      iwbeg = nint((wbeg-b)/delta)+1
      iwend = nint((wend-b)/delta)+1
      iwbeg=f*(iwbeg-1) + 1
      iwend=f*(iwend-1) + 1 + f*np2*itlag_step
      call zwindow(xinterp,yinterp,ninterp,np,iwbeg,iwend,xw,yw,nw)
      do i=1,nw
         timew(i) = wbeg + delta*real(i-1)/real(f)
      enddo
c*****CALC UNCORRECTED FAST AND SLOW COMPONENTS **
      call zrotate2d(xw,yw,nw,np,fast,ssw,sfw)
c*****CALC CORRECTED FAST AND SLOW COMPONENT**
      call zlag(ssw,sfw,nw,np,ilag,np2int,sscw,sfcw,noverlap)

c  ** WRITE OUT PARTICLE MOTION **
c  ** rotate corrected fast and slow waves back to original coords **
      call zrotate2d(sscw,sfcw,noverlap,np,-fast,xcw,ycw)
c  ** the points we required are numbered 1-noverlap **
c  ** find normalising factor for the particle motion and normalise to
c      a maximum value of 1**
      call zabs_max(xw,noverlap,np,max_xw)
      call zabs_max(yw,noverlap,np,max_yw)
      call zabs_max(xcw,noverlap,np,max_xcw)
      call zabs_max(ycw,noverlap,np,max_ycw)
      norm_pm=1./max(max_xw,max_yw,max_xcw,max_ycw)
      call zmultiply(xw,noverlap,np,norm_pm,xwn)
      call zmultiply(yw,noverlap,np,norm_pm,ywn)
      call zmultiply(xcw,noverlap,np,norm_pm,xcwn)
      call zmultiply(ycw,noverlap,np,norm_pm,ycwn)
c*****UNCORRECTED PARTICLE MOTION IN S-WAVE WINDOW **
c  ** particle motion is composed of xwn,ywn in window wbeg-wend **
      call zwrite2(file_pm,lu,fmt2,xwn,ywn,noverlap,np)
c*****CORRECTED PARTICLE MOTION IN S-WAVE WINDOW **
c  ** particle motion is composed of xcwnw,ycwn in window wbeg-wend **
      call zwrite2(file_pmc,lu,fmt2,xcwn,ycwn,noverlap,np)

c  ** WRITE OUT FAST AND SLOW WAVEFORMS **
c  ** normalise the S-waves **
      call zabs_max(ssw,noverlap,np,max_ssw)
      call zabs_max(sfw,noverlap,np,max_sfw)
      call zabs_max(sscw,noverlap,np,max_sscw)
      call zabs_max(sfcw,noverlap,np,max_sfcw)
      call zmultiply(sfw,noverlap,np,1./max_sfw,sfwn)
      call zmultiply(sfcw,noverlap,np,1./max_sfcw,sfcwn)
      call zcovariance(sscw,sfcw,noverlap,np,cov)
c  ** if cross correlation is negative then flip the slow wave so that
c      the match of the waveforms is more obvious **
      if (cov(1,2).ge.0.) then
         call zmultiply(ssw,noverlap,np,1./max_ssw,sswn)
         call zmultiply(sscw,noverlap,np,1./max_sscw,sscwn)
      else
         call zmultiply(ssw,noverlap,np,-1./max_ssw,sswn)
         call zmultiply(sscw,noverlap,np,-1./max_sscw,sscwn)
      endif
c*****UNCORRECTED FAST AND SLOW COMPONENTS **
      call zwrite2(file_ss,lu,fmt2,timew,sswn,noverlap,np)
      call zwrite2(file_sf,lu,fmt2,timew,sfwn,noverlap,np)
c*****CORRECTED FAST AND SLOW COMPONENTS **
      call zwrite2(file_ssc,lu,fmt2,timew,sscwn,noverlap,np)
      call zwrite2(file_sfc,lu,fmt2,timew,sfcwn,noverlap,np)

c*****INTERPOLATED ERROR SURFACE **
      open(lu,file=file_error)
      do i=1,np1
         do j=1,np2int
           test_tlag = delta*real((j-1)*itlag_step)/real(f)
           test_fast = -90. + 180.*real(i-1)/real(np1-1)
           write(lu,*) test_tlag,test_fast,error_int(i,j)
         enddo
      enddo
      close(lu)

c--------- PLOTTING PARAMETERS REQUIRED BY GMT --------
      open(lu,file=file_gmt)
      write(lu,*) 'ppick ',ppick
      write(lu,*) 'spick ',spick
      write(lu,*) 'wbeg ',wbeg
      write(lu,*) 'wend ',wend
      write(lu,*) 'waveform_minortick ',ftick((wend-wbeg)/4.)
      write(lu,*) 'waveform_majortick ',ftick((wend-wbeg))
      write(lu,*) 'tlag_scale ',tlag_scale
      write(lu,*) 'tlag_minortick ',ftick(tlag_scale/4.)
      write(lu,*) 'tlag_majortick ',ftick(tlag_scale)
      write(lu,*) 'fast_scale ',fast_scale
      write(lu,*) 'delta ',delta
      write(lu,*) 'error_grid_tlag_int ',itlag_step*delta/real(f)
      write(lu,*) 'rt_max_y',max_rt
      write(lu,*) 'rt_minortick_y',ftick(max_rt/4.)
      write(lu,*) 'rt_majortick_y',ftick(max_rt)
C     * set up gmt trace plotting window
      gmt_trace_start = wbeg - 0.2
      if (gmt_trace_start.lt.b) gmt_trace_start = b
      gmt_trace_end = wend + 0.2
      if (gmt_trace_end.gt.(real(n-1)*delta))
     +     gmt_trace_end = real(n-1)*delta

      write(lu,*) 'rt_b_x',gmt_trace_start
      write(lu,*) 'rt_e_x',gmt_trace_end
      write(lu,*) 'rt_minortick_x',ftick((real(n-1)*delta)/32.)
      write(lu,*) 'rt_majortick_x',ftick((real(n-1)*delta)/8.)
      close(lu)

      return
      end

c-----------------------------------------------------------------------
      function ftick(range)
c-----------------------------------------------------------------------
      implicit none
      real range,ftick
      if (range.le.0.0078125) then
        ftick=0.001
      else if (range.le.0.015625) then
        ftick=0.0025
      else if (range.le.0.03125) then
        ftick=0.005
      else if (range.le.0.0625) then
        ftick=0.01
      else if (range.le.0.125) then
        ftick=0.025
      else if (range.le.0.25) then
        ftick=0.05
      else if (range.le.0.5) then
        ftick=0.1
      else if (range.le.1.0) then
        ftick=0.2
      else if (range.le.2.0) then
        ftick=0.5
      else if (range.le.4.0) then
        ftick=1.
      else if (range.le.8.0) then
        ftick=2.
      else if (range.le.16.0) then
        ftick=5.
      else if (range.le.32.0) then
        ftick=10.
      else if (range.le.64.0) then
        ftick=20.
      else if (range.le.128.0) then
        ftick=50.
      else if (range.le.256.0) then
        ftick=50.
      else if (range.le.512.0) then
        ftick=100.
      else if (range.le.1024.0) then
        ftick=200.
      else if (range.le.2048.0) then
        ftick=500.
      else if (range.le.4096.0) then
        ftick=1000.
      else if (range.le.8192.0) then
        ftick=2000.
      else if (range.le.16384.0) then
        ftick=5000.
      else if (range.le.32768.0) then
        ftick=10000.
      else
        ftick=range
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine zwritesac(file,lu,np,h1,h2,h3a,h3b,h4,d1,npts)
c-----------------------------------------------------------------------
c
c      write a sac ascii data file
c        write header info (in 4 parts)
c        write data part1
c        (may need to add a data part2 if unevenly spaced data is ever used)
c
c      variables:
c input: file      char*50            file to write output to
c         lu            int                  logical unit to read file to
c         np            int                  dimension of d1 array
c         h1            real(14,5)            header part1
c         h2            int(8,5)            header part2
c         h3a      char*8            header part3a
c         h3b      char*16            header part3b
c         h4            char*8 (7,3)      header part4
c         d1            real(np)            data points
c         npts      int                  number of data points (check with h2(2,5))
c
c-----------------------------------------------------------------------
c  modifications:
c      11-07-01      N. Teanby      Original code
c-----------------------------------------------------------------------

      implicit none
      integer lu,np
      real h1(14,5)
      integer h2(8,5)
      character h3a*8, h3b*16
      character*8 h4(7,3)
      real d1(np)
      integer i,npts
      character*50 file

c  ** open data file **
      open(lu,file=file,status='unknown')

c  ** write header info **
      do 1 i=1,14
         write(lu,1001) h1(i,1),h1(i,2),h1(i,3),h1(i,4),h1(i,5)
1      continue
      do 2 i=1,8
         write(lu,1002) h2(i,1),h2(i,2),h2(i,3),h2(i,4),h2(i,5)
2      continue
      write(lu,1003) h3a,h3b
      do 4 i=1,7
         write(lu,1004) h4(i,1),h4(i,2),h4(i,3)
4      continue

c  ** check npts with npts from header **
      if (npts.ne.h2(2,5)) then
            print*,'WARNING:zsacwrite'
            print*,'WARNING:  npts supplied .ne. npts from header'
      endif
c  ** check npts does not exceed array dimension **
      if (npts.gt.np) then
         print*,'WARNING:zsacwrite'
         print*,'WARNING:  npts.gt.np, data missed out'
      endif

c  ** remove v small numbers or SAC complains !!! **
      do 50 i=1,npts
         if (abs(d1(i)).le.1.e-14) then
            d1(i)=0.
         endif
50      continue

c  ** write data **
      do 100 i=1,5*int(npts/5),5
         write(lu,1010) d1(i),d1(i+1),d1(i+2),d1(i+3),d1(i+4)
100      continue
c  ** check that npts is a multiple of 5 **
      if (int(npts/5).ne.npts/5) then
         print*,'WARNING:zsacwrite'
         print*,'WARNING:  npts not multiple of 5'
         print*,'WARNING:  last few (<5) points may be missed off'
      endif

      close(lu)
      return

c **      format statements **
1001      format(5g15.7)
1002      format(5i10)
1003      format(a8,a16)
1004      format(3a8)
1010      format(5g15.7)

      end

      function fstrcat(str1,str2)
c      concatenate two strings, str1 and str2
      implicit none
      character*50 str1,str2,fstrcat
      integer len1,len2

      len1=index(str1,' ')-1
      len2=index(str2,' ')-1
      if (len1+len2.gt.50) then
         print*,'string1 = ',str1
         print*,'string2 = ',str2
         pause 'ERROR: fstrcat: combined string length exceeds 50'
      else
         fstrcat=str1(1:len1)//str2(1:len2)
      endif
      return
      end

c-----------------------------------------------------------------------
      function fftable(ndf)
c-----------------------------------------------------------------------
c
c      function to return value of f statistic
c
c      specialised for k=2 and alpha=0.05 (2 params and 95% confidence)
c
c      table downloaded from:
c      http://www.itl.nist.gov/div898/handbook/eda/section3/eda3673.htm
c
c      ndf       int      number of degrees of freedom
c
c-----------------------------------------------------------------------
c      N. Teanby      1-8-02      Original code
c-----------------------------------------------------------------------

      implicit none
      real fftable
      integer ndf
      real ftable_data(100)
      data ftable_data / 199.500, 19.000, 9.552, 6.944, 5.786,
     >5.143, 4.737, 4.459, 4.256, 4.103,
     >3.982, 3.885, 3.806, 3.739, 3.682,
     >3.634, 3.592, 3.555, 3.522, 3.493,
     >3.467, 3.443, 3.422, 3.403, 3.385,
     >3.369, 3.354, 3.340, 3.328, 3.316,
     >3.305, 3.295, 3.285, 3.276, 3.267,
     >3.259, 3.252, 3.245, 3.238, 3.232,
     >3.226, 3.220, 3.214, 3.209, 3.204,
     >3.200, 3.195, 3.191, 3.187, 3.183,
     >3.179, 3.175, 3.172, 3.168, 3.165,
     >3.162, 3.159, 3.156, 3.153, 3.150,
     >3.148, 3.145, 3.143, 3.140, 3.138,
     >3.136, 3.134, 3.132, 3.130, 3.128,
     >3.126, 3.124, 3.122, 3.120, 3.119,
     >3.117, 3.115, 3.114, 3.112, 3.111,
     >3.109, 3.108, 3.107, 3.105, 3.104,
     >3.103, 3.101, 3.100, 3.099, 3.098,
     >3.097, 3.095, 3.094, 3.093, 3.092,
     >3.091, 3.090, 3.089, 3.088, 3.087 /

      if (ndf.le.0) then
c      ** if ndf is below range 1-100 error **
         pause 'ERROR: fftable: ndf.le.0'
      else if (ndf.le.100) then
c      ** if ndf is in range 1-100 take directly from tabulated values **
         fftable=ftable_data(ndf)
      else if (ndf.le.999) then
c      ** if ndf is in range 100-999 use linear interpolation **
c      ** f=3.087 at ndf = 100, and f=3.000 at ndf=999 **
         fftable = 3.087 - 0.087*real(ndf-100)/899.
      else if (ndf.gt.999) then
c      ** if ndf is over range return value at 999 **
         fftable = 3.0
         print*,'WARNING: ndf.gt.999, f statistic at ndf=999 returned'
      endif

      return
      end

c-----------------------------------------------------------------------
      function fda(angle1,angle2)
c-----------------------------------------------------------------------
c      find differnce between two angles which have 180 deg ambiguity
c      angle1/2      angles in degrees
c      fda            difference in degrees
c-----------------------------------------------------------------------
c      N. Teanby      30-8-02      Original code
c-----------------------------------------------------------------------
      implicit none
      real angle1,angle2,pi,fda
      parameter (pi=3.141592654)

      fda = acos( abs(cos( (angle1-angle2)*pi/180.) ) )*180./pi

      return
      end

