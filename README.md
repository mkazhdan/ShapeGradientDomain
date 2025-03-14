<CENTER><H2>Shape Gradient Domain (Version 5.0)</H2></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLE">executable</A>
<A HREF="#NOTES">notes</A>
<A HREF="#EXAMPLES">examples</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
[<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">Chuang et al. 2009</A>], [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">Prada et al. 2016</A>], [<A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">Chuang et al. 2016</A>]<br>
<A HREF="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.x64.zip">Windows (x64) Executables</A><BR>
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.zip">Source Code</A><br>
<A HREF="https://github.com/mkazhdan/ShapeGradientDomain">GitHub Repository</A><BR>
(Older Versions:
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.5/">V4.5</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.0/">V4.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version3.0/">V3.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version2.0/">V2.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version1.0/">V1.0</A>)
<br>


<HR>
<A NAME="DESCRIPTION"><B>DESCRIPTION</B></A><br>
<UL>
<P>
<!--<LI> <B><U>ShapeGradientDomain</U></B>-->

<LI>This code performs gradient domain processing [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">Chuang et al. 2009</A>, <A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">Chuang et al. 2016</A>] on signals defined on a mesh, where the signal can be either a color-field represented as a color per vertex or is the position of the vertices themselves.
The code supports both sharpening and smoothing of the signals through the solution of a screened-Poisson equation.
Specifically, given an input signal <I>F</I>, it solves for the signal <I>G</I> minimizing:<BR>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &beta;&sdot;||&lambda;&sdot;&nabla;<I>F</I> - &nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &beta; is the gradient-fitting weight and &lambda; is the gradient scale factor.

<LI> The code supports inhomogenous processing by allowing the user to replace the Riemannian metric, <I>g</I>, given by the embedding, with a metric that adjusts to the curvature of the surface. Specifically, given orthonormal principal curvature directions, the (idenity) metric is replaced with:<BR>
<CENTER> <i>g</I> &larr; Id. + &epsilon;&sdot;&Kappa;<sup>2</sup></CENTER>
where Id. is the identity matrix, &Kappa;<sup>2</sup> is the diagonal matrix whose entries are the sum of the squares of the principal curature values, and &epsilon; is the curvature weight. (When &epsilon;=0, this reproduces the Riemannian metric given by the embedding.)<br>
Curvatures are estimated using the surface normals. If none are provided, the vertex normals are estimated as the area-weighted sum of adjacent triangle normals.

<LI> For robustness, the implementation supports pre-smoothing the normals before estimating curvatures. Similar to above, this amounts to minimizing an energy of the form:<br>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &gamma;&sdot;||&nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &gamma; is the normal smoothing weight.<br>

<LI>In the case that the signal processed describes the vertices' normals, we also support fitting the original geometry to the processed normals using the approach of [<A HREF="https://dl.acm.org/doi/10.1145/1015706.1015774">Yu et al. 2004</A>].<BR>
 Given input positions <I>F</I> and a target normal field <I>N</I> This amounts to minimizing:<br>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &delta;&sdot;||(&nabla;<I>F</I> - <I>N</I>&sdot;&#10216;&nabla;<I>F</I>,<I>N</I>&#10217;) - &nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &delta; is the projected gradient fitting weight weight.<br>

</UL>
<!--
Note that when &beta; is set to zero the two executables differ in that the first emaulates harmonic flow from the input geometry to Euclidean three-space (allowing the variation at a vertex to occur in any direction) while the second emulates harmonic flow from the input geometry to the two-sphere (constrainting the variation at a vertex to occur in the tangent space of the associated normal).
-->

<HR>
<a name="EXECUTABLE"><b>EXECUTABLE</b></a><br>


<UL>
<DL>
<DETAILS>
<SUMMARY>
<font size="+1"><b>ShapeGradientDomain</b></font>:
Processes either the vertex positions, or per-vertex colors, performing homogeneous/inhomogeneous gradient-domain smoothing and sharpening (with the option of smoothing normals in a pre-processing step)
</SUMMARY>

<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--nWeight</b> &#60;<i>normal smoothing time</i>&#62;]
<DD> This floating point value gives the normal smoothing weight (&gamma;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--gWeight</b> &#60;<i>gradient interpolation weight</i>&#62;]
<DD> This floating point value gives the weight for gradient interpolation (&beta;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--gScale</b> &#60;<i>gradient scale</i>&#62;]
<DD> This floating point value gives the scale factor for the target gradient field (&lambda;).<BR>
The default value for this parameter is 1.0.<br>

<DT>[<b>--kWeight</b> &#60;<i>curvature weight</i>&#62;]
<DD> This floating point value gives the curvature weight for adjusting the metric (&epsilon;).<BR>
The default value for this parameter is 0.0.<br>

<DT>[<b>--npWeight</b> &#60;<i>normal projection weight</i>&#62;]
<DD> This floating point value gives the normal projection weight (&delta;).<BR>
The default value for this parameter is 10<SUP>2</SUP>.<br>

<DT>[<b>--signal</b>]
<DD> This inter value describes which signal is to be processed:
<UL>
<LI><B>0</B>: Vertex positions
<LI><B>1</B>: Vertex normals (will be assigned if none are give)
<LI><B>2</B>: Vertex colors
<LI><B>3</B>: Vertex normals used to fit the geometry.
</UL>
The default value for this parameter is 0.

<DT>[<b>--verbose</b>]
<DD> If this flag is enabled, the code will output processing information.

</DETAILS>

<DETAILS>
<SUMMARY>
<font size="+1"><b>ShapeCurvature</b></font>:
Visualizes the estimated total curvature as a gray-scale per-vertex color.
</SUMMARY>

<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--nWeight</b> &#60;<i>normal smoothing weight</i>&#62;]
<DD> This floating point value gives the normal smoothing weight (&gamma;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--verbose</b>]
<DD> If this flag is enabled, the code will output processing information.

</DETAILS>
</DL>
</UL>



<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>
<UL>
<LI> The code requires <A HREF="https://eigen.tuxfamily.org">Eigen</A> as a numerical solver.<BR>
If you are using Eigen and your implementation is backed by <A HREF="https://software.intel.com/en-us/intel-mkl/">Intel's Math Kernel Library</A> (see discussion <A HREF="https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html">here</A>), enable the <CODE>EIGEN_USE_MKL_ALL</CODE> macro by defining it in the file <CODE>PreProcessor.h</CODE>. (The two versions of the <A HREF="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.x64.zip">Windows executables</A> are compiled with and without MKL support.)
<LI> For consistency of parameters, surfaces are rescaled to have unit area prior to processing (and then scaled back after processing).
</UL>

<HR>
<A NAME="EXAMPLES"><B>EXAMPLES</B></A><br>
<UL>
<LI><B>ShapeGradientDomain</B>:
<P>The figure below shows example of both homogenous and inhomogeneous geometry processing.<BR>
<UL>
<LI><B>Top row</B>: Homogeneous sharpening:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.sharp.ply --gScale 2 --gWeight &lt;gradient weight&gt; </CODE></PRE>
<LI><B>Middle row</B>: Homogeneous smoothing:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.smooth.ply --gScale 0 --gWeight &lt;gradient weight&gt; </CODE></PRE>
<LI><B>Bottom row</B>: Inhomogeneous smoothing:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.smooth.ply --gScale 0 --gWeight &lt;gradient weight&gt; --kWeight 0.02</CODE></PRE>
</UL>
For all three, the value of <CODE>&lt;gradient weight&gt;</CODE> is chosen from {<CODE>1e-3</CODE>, <CODE>1e-4</CODE>, <CODE>1e-5</CODE>}
<CENTER><IMG SRC="armadillo.png" WIDTH=80%></CENTER>

<LI><B>ShapeCurvature</B>:
<P>The figure below visualizes the total curvature (square root of the sum of the squares of the principal curvatures) as per-vertex grayscale colors and gives the standard deviations of the estimated total curvature &sigma;<SUB>|&kappa;|</SUB> 
<PRE>     <CODE>ShapeCurvature --in armadillo.ply --out armadillo.curvature.ply --nWeight &lt;normal smoothing weight&gt;</CODE></PRE>
The value of <CODE>&lt;normal smoothing weight&gt;</CODE> is chosen from {<CODE>1e-4</CODE>,<CODE>1e-5</CODE>,<CODE>1e-6</CODE>}.<BR>
<CENTER><IMG SRC="armadillo.curvature.png" WIDTH=80%></CENTER>
</UL>




<HR>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version2.0/">Version 2.0</a>:
<OL>
<LI> Added options to weight both value and gradient interpolation terms.
<LI> Changed default values to correspond to diffusion time.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version3.0/">Version 3.0</a>:
<OL>
<LI> Integrated normal smoothing within the gradient-domain processing application.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.0/">Version 4.0</a>:
<OL>
<LI> Fixed default value weight to 1.0.
<LI> Added support for fitting the geometry to prescribed normals.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.5/">Version 4.5</a>:
<OL>
<LI> Enabled symbolic matrix factorization in pre-processing.
<LI> Added support (and made default) normal smoothing through per-channel diffusion.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/">Version 5.0</a>:
<OL>
<LI> Removed harmonic normal smoothing
<LI> Removed normal smoothing iterations
<LI> Added cuvature visualization executable
</OL>

</DETAILS><CENTER><H2>Shape Gradient Domain (Version 5.0)</H2></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLE">executable</A>
<A HREF="#NOTES">notes</A>
<A HREF="#EXAMPLES">examples</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
[<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">Chuang et al. 2009</A>], [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">Prada et al. 2016</A>], [<A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">Chuang et al. 2016</A>]<br>
<A HREF="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.x64.zip">Windows (x64) Executables</A><BR>
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.zip">Source Code</A><br>
<A HREF="https://github.com/mkazhdan/ShapeGradientDomain">GitHub Repository</A><BR>
(Older Versions:
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.5/">V4.5</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.0/">V4.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version3.0/">V3.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version2.0/">V2.0</A>,
<A href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version1.0/">V1.0</A>)
<br>


<HR>
<A NAME="DESCRIPTION"><B>DESCRIPTION</B></A><br>
<UL>
<P>
<!--<LI> <B><U>ShapeGradientDomain</U></B>-->

<LI>This code performs gradient domain processing [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">Chuang et al. 2009</A>, <A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">Chuang et al. 2016</A>] on signals defined on a mesh, where the signal can be either a color-field represented as a color per vertex or is the position of the vertices themselves.
The code supports both sharpening and smoothing of the signals through the solution of a screened-Poisson equation.
Specifically, given an input signal <I>F</I>, it solves for the signal <I>G</I> minimizing:<BR>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &beta;&sdot;||&lambda;&sdot;&nabla;<I>F</I> - &nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &beta; is the gradient-fitting weight and &lambda; is the gradient scale factor.

<LI> The code supports inhomogenous processing by allowing the user to replace the Riemannian metric, <I>g</I>, given by the embedding, with a metric that adjusts to the curvature of the surface. Specifically, given orthonormal principal curvature directions, the (idenity) metric is replaced with:<BR>
<CENTER> <i>g</I> &larr; Id. + &epsilon;&sdot;&Kappa;<sup>2</sup></CENTER>
where Id. is the identity matrix, &Kappa;<sup>2</sup> is the diagonal matrix whose entries are the sum of the squares of the principal curature values, and &epsilon; is the curvature weight. (When &epsilon;=0, this reproduces the Riemannian metric given by the embedding.)<br>
Curvatures are estimated using the surface normals. If none are provided, the vertex normals are estimated as the area-weighted sum of adjacent triangle normals.

<LI> For robustness, the implementation supports pre-smoothing the normals before estimating curvatures. Similar to above, this amounts to minimizing an energy of the form:<br>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &gamma;&sdot;||&nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &gamma; is the normal smoothing weight.<br>

<LI>In the case that the signal processed describes the vertices' normals, we also support fitting the original geometry to the processed normals using the approach of [<A HREF="https://dl.acm.org/doi/10.1145/1015706.1015774">Yu et al. 2004</A>].<BR>
 Given input positions <I>F</I> and a target normal field <I>N</I> This amounts to minimizing:<br>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &delta;&sdot;||(&nabla;<I>F</I> - <I>N</I>&sdot;&#10216;&nabla;<I>F</I>,<I>N</I>&#10217;) - &nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &delta; is the projected gradient fitting weight weight.<br>

</UL>
<!--
Note that when &beta; is set to zero the two executables differ in that the first emaulates harmonic flow from the input geometry to Euclidean three-space (allowing the variation at a vertex to occur in any direction) while the second emulates harmonic flow from the input geometry to the two-sphere (constrainting the variation at a vertex to occur in the tangent space of the associated normal).
-->

<HR>
<a name="EXECUTABLE"><b>EXECUTABLE</b></a><br>


<UL>
<DL>
<DETAILS>
<SUMMARY>
<font size="+1"><b>ShapeGradientDomain</b></font>:
Processes either the vertex positions, or per-vertex colors, performing homogeneous/inhomogeneous gradient-domain smoothing and sharpening (with the option of smoothing normals in a pre-processing step)
</SUMMARY>

<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--nWeight</b> &#60;<i>normal smoothing time</i>&#62;]
<DD> This floating point value gives the normal smoothing weight (&gamma;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--gWeight</b> &#60;<i>gradient interpolation weight</i>&#62;]
<DD> This floating point value gives the weight for gradient interpolation (&beta;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--gScale</b> &#60;<i>gradient scale</i>&#62;]
<DD> This floating point value gives the scale factor for the target gradient field (&lambda;).<BR>
The default value for this parameter is 1.0.<br>

<DT>[<b>--kWeight</b> &#60;<i>curvature weight</i>&#62;]
<DD> This floating point value gives the curvature weight for adjusting the metric (&epsilon;).<BR>
The default value for this parameter is 0.0.<br>

<DT>[<b>--npWeight</b> &#60;<i>normal projection weight</i>&#62;]
<DD> This floating point value gives the normal projection weight (&delta;).<BR>
The default value for this parameter is 10<SUP>2</SUP>.<br>

<DT>[<b>--signal</b>]
<DD> This inter value describes which signal is to be processed:
<UL>
<LI><B>0</B>: Vertex positions
<LI><B>1</B>: Vertex normals (will be assigned if none are give)
<LI><B>2</B>: Vertex colors
<LI><B>3</B>: Vertex normals used to fit the geometry.
</UL>
The default value for this parameter is 0.

<DT>[<b>--verbose</b>]
<DD> If this flag is enabled, the code will output processing information.

</DETAILS>

<DETAILS>
<SUMMARY>
<font size="+1"><b>ShapeCurvature</b></font>:
Visualizes the estimated total curvature as a gray-scale per-vertex color.
</SUMMARY>

<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--nWeight</b> &#60;<i>normal smoothing weight</i>&#62;]
<DD> This floating point value gives the normal smoothing weight (&gamma;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--verbose</b>]
<DD> If this flag is enabled, the code will output processing information.

</DETAILS>
</DL>
</UL>



<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>
<UL>
<LI> The code requires <A HREF="https://eigen.tuxfamily.org">Eigen</A> as a numerical solver.<BR>
If you are using Eigen and your implementation is backed by <A HREF="https://software.intel.com/en-us/intel-mkl/">Intel's Math Kernel Library</A> (see discussion <A HREF="https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html">here</A>), enable the <CODE>EIGEN_USE_MKL_ALL</CODE> macro by defining it in the file <CODE>PreProcessor.h</CODE>. (The two versions of the <A HREF="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/ShapeGradientDomain.x64.zip">Windows executables</A> are compiled with and without MKL support.)
<LI> For consistency of parameters, surfaces are rescaled to have unit area prior to processing (and then scaled back after processing).
</UL>

<HR>
<A NAME="EXAMPLES"><B>EXAMPLES</B></A><br>
<UL>
<LI><B>ShapeGradientDomain</B>:
<P>The figure below shows example of both homogenous and inhomogeneous geometry processing.<BR>
<UL>
<LI><B>Top row</B>: Homogeneous sharpening:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.sharp.ply --gScale 2 --gWeight &lt;gradient weight&gt; </CODE></PRE>
<LI><B>Middle row</B>: Homogeneous smoothing:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.smooth.ply --gScale 0 --gWeight &lt;gradient weight&gt; </CODE></PRE>
<LI><B>Bottom row</B>: Inhomogeneous smoothing:
<PRE>     <CODE>ShapeGradientDomain --in armadillo.ply --out armadillo.smooth.ply --gScale 0 --gWeight &lt;gradient weight&gt; --kWeight 0.02</CODE></PRE>
</UL>
For all three, the value of <CODE>&lt;gradient weight&gt;</CODE> is chosen from {<CODE>1e-3</CODE>, <CODE>1e-4</CODE>, <CODE>1e-5</CODE>}
<CENTER><IMG SRC="armadillo.png" WIDTH=80%></CENTER>

<LI><B>ShapeCurvature</B>:
<P>The figure below visualizes the total curvature (square root of the sum of the squares of the principal curvatures) as per-vertex grayscale colors and gives the standard deviations of the estimated total curvature &sigma;<SUB>|&kappa;|</SUB> 
<PRE>     <CODE>ShapeCurvature --in armadillo.ply --out armadillo.curvature.ply --nWeight &lt;normal smoothing weight&gt;</CODE></PRE>
The value of <CODE>&lt;normal smoothing weight&gt;</CODE> is chosen from {<CODE>1e-4</CODE>,<CODE>1e-5</CODE>,<CODE>1e-6</CODE>}.<BR>
<CENTER><IMG SRC="armadillo.curvature.png" WIDTH=80%></CENTER>
</UL>




<HR>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version2.0/">Version 2.0</a>:
<OL>
<LI> Added options to weight both value and gradient interpolation terms.
<LI> Changed default values to correspond to diffusion time.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version3.0/">Version 3.0</a>:
<OL>
<LI> Integrated normal smoothing within the gradient-domain processing application.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.0/">Version 4.0</a>:
<OL>
<LI> Fixed default value weight to 1.0.
<LI> Added support for fitting the geometry to prescribed normals.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version4.5/">Version 4.5</a>:
<OL>
<LI> Enabled symbolic matrix factorization in pre-processing.
<LI> Added support (and made default) normal smoothing through per-channel diffusion.
</OL>

<a href="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version5.0/">Version 5.0</a>:
<OL>
<LI> Removed harmonic normal smoothing
<LI> Removed normal smoothing iterations
<LI> Added cuvature visualization executable
</OL>

</DETAILS>