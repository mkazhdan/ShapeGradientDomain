<CENTER><H2>Shape Gradient Domain (Version 2.0)</H2></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLES">executables</A>
<A HREF="#EXAMPLES">examples</A>
<A HREF="#NOTES">notes</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">SGP 2009 Paper</A>, <A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">SIGGRAPH 2016 Paper</A>, <A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">JCGT 2016 Paper</A><br>
<A HREF="https://www.cs.jhu.edu/~misha/ShapeGradientDomain/Version2.0/ShapeGradientDomain.x64.zip">Windows (x64) Executables</A><BR>
<A href="https://www.cs.jhu.edu/~misha/ShapeGradientDomain/Version2.0/ShapeGradientDomain.zip">Source Code</A><br> <A HREF="https://github.com/mkazhdan/ShapeGradientDomain">GitHub Repository</A><BR>
(Older Versions:
<A href="https://www.cs.jhu.edu/~misha/ShapeGradientDomain/Version1.0/">V1.0</A>)
<br>


<HR>
<A NAME="DESCRIPTION"><B>DESCRIPTION</B></A><br>
The code is comprised of two executables.
<UL>
<P>
<LI> <B><U>ShapeGradientDomain</U></B>: [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">SGP 2009</A>]<P>
This code performs gradient domain processing on signals defined on a mesh, where the signal can be either a color-field represented as a color per vertex or is the position of the vertices themselves..
The code supports both sharpening and smoothing of the signals through the solution of a screened-Poisson equation.
Specifically, given an input signal <I>F</I>, it solves for the signal <I>G</I> minimizing:<BR>
<CENTER>
<I>E</I>(<I>G</I>) = &alpha;&sdot;||<I>F</I>-<I>G</I>||<sup>2</sup> + &beta;&sdot;||&lambda;&sdot;&nabla;<I>F</I> - &nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &alpha; is the value-fitting weight, &beta; is the gradient-fitting weight, and &lambda; is the gradient scale factor.<br>
The code supports inhomogenous processing by allowing the user to replace the Riemannian metric, <I>g</I>, given by the embedding, with a metric that adjusts to the curvature of the surface. Specifically, given orthonormal principal curvature directions, the (idenity) metric is replaced with:<BR>
<CENTER> Id. + &epsilon;&sdot;&Kappa;<sup>2</sup></CENTER>
where Id. is the identity matrix and &Kappa;<sup>2</sup> is the diagonal matrix whose entries are the squares of the principal curature values and &epsilon; is the curvature weight.<br>
Curvatures are estimated using the surface normals. If none are provided, the vertex normals are estimated as the area-weighted sum of adjacent triangle normals.
<LI> <B><U>Normal Smooth</U></B>: [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">SIGGRAPH 2016</A>]<P>
This code performs multiple iterations of harmonic smoothing of the surface normals. As with the code above, this amounts to minimizing:<br>
<CENTER>
<I>E</I>(<I>G</I>) = ||<I>F</I>-<I>G</I>||<sup>2</sup> + &gamma;&sdot;||&nabla;<I>G</I>||<sup>2</sup>
</CENTER><BR>
where &gamma; is the diffusion weight (time).<br>

If no normals are provided, the vertex normals are estimated as the area-weighted sum of adjacent triangle normals.
</UL>
Note that when &beta; is set to zero the two executables differ in that the first emaulates harmonic flow from the input geometry to Euclidean three-space (allowing the variation at a vertex to occur in any direction) while the second emulates harmonic flow from the input geometry to the two-sphere (constrainting the variation at a vertex to occur in the tangent space of the associated normal).

<HR>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>


<UL>
<DL>
<DETAILS>
<SUMMARY>
<font size="+1"><b>ShapeGradientDomain</b></font>:
Processes either the vertex positions, or per-vertex colors, performing isotropic/anisotropic gradient-domain smoothing and sharpening
[<A href="https://www.cs.jhu.edu/~misha/MyPapers/SGP09.pdf">SGP 2009</A>, <A href="https://www.cs.jhu.edu/~misha/MyPapers/JCGT16.pdf">JCGT 2016</A>]
</SUMMARY>

<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--vWeight</b> &#60;<i>value interpolation weight</i>&#62;]
<DD> This floating point value gives the weight for value interpolation (&alpha;).<BR>
The default value for this parameter is 1.<br>

<DT>[<b>--gWeight</b> &#60;<i>gradient interpolation weight</i>&#62;]
<DD> This floating point value gives the weight for gradient interpolation (&beta;).<BR>
The default value for this parameter is 10<sup>-4</sup>.<br>

<DT>[<b>--gScale</b> &#60;<i>gradient scale</i>&#62;]
<DD> This floating point value gives the scale factor for the target gradient field (&lambda;).<BR>
The default value for this parameter is 1.0.<br>

<DT>[<b>--kWeight</b> &#60;<i>curvature weight</i>&#62;]
<DD> This floating point value gives the curvature weight for adjusting the metric (&epsilon;).<BR>
The default value for this parameter is 0.0.<br>

<DT>[<b>--useColors</b>]
<DD> If this flag is enabled, the signal to be processed is the per-vertex color field. Otherwise, it is the vertex positions.<BR>
If the flag is enabled and the input file does not contain per-vertex colors, colors will be assigned from the normals.

<DT>[<b>--verbose</b>]
<DD> If this flag is enabled, the code will output processing information.

</DETAILS>
</DL>
</UL>



<UL>
<DL>
<DETAILS>
<SUMMARY>
<font size="+1"><b>NormalSmooth</b></font>:
Diffuses surface normals, restricting the change to be within the tangent plane [<A href="https://www.cs.jhu.edu/~misha/MyPapers/SIG16.pdf">SIGGRAPH 2016</A>].
</SUMMARY>


<DT><b>--in</b> &#60;<i>input geometry</i>&#62;
<DD> This string specifies the name of the input geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>ouput geometry</i>&#62;]
<DD> This string specifies the name of the output geometry, represented in <A HREF="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--iters</b> &#60;<i>iterations</i>&#62;]
<DD> This integer value specifies the number of smoothing iterations that are to be performed..<BR>
The default value for this parameter is 1.<br>

<DT>[<b>--dTime</b> &#60;<i>gradient interpolation weight / diffusion time</i>&#62;]
<DD> This floating point value gives the weight for gradient interpolation / diffusion time (&gamma;).<BR>
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
If you are using Eigen and your implementation is backed by <A HREF="https://software.intel.com/en-us/intel-mkl/">Intel's Math Kernel Library</A> (see discussion <A HREF="https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html">here</A>), enable the <CODE>EIGEN_USE_MKL_ALL</CODE> macro by defining it in the file <CODE>PreProcessor.h</CODE>. (The two versions of the <A HREF="https://www.cs.jhu.edu/~misha/Code/ShapeGradientDomain/Version2.0/ShapeGradientDomain.x64.zip">Windows executables</A> are similarly compiled with and without MKL support.)
</UL>

<HR>
<A NAME="EXAMPLES"><B>EXAMPLES</B></A><br>
The figure below shows example of both isotropic and anisotropic geometry processing.<BR>
<UL>
<LI> Geometric effects are obtained by either amplifying (&lambda;=2 in the top row) or dampening (&lambda;=0 in the bottom two rows) gradients.
<LI> From left to right, the gradient weight is successively decreased (&beta;=10<sup>-3</sup>, &beta;=10<sup>-4</sup>, and &beta;=10<sup>-5</sup>) corresponding to successively more loacalized edits.
<LI> The processing is isotropic in the top two rows (&epsilon;=0 and &gamma; is irrelevent) and isotropic in the bottom one (&epsilon;=0.02 and &gamma;=10<sup>-4</sup>).
</UL>

<CENTER><IMG SRC="armadillo.png" WIDTH=80%></CENTER>



<HR>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>
<a href="../Version2.0/">Version 2.0</a>:
<OL>
<LI> Added options to weight both value and gradient interpolation terms.
<LI> Changed default values to correspond to diffusion time.
</OL>

</DETAILS>
