/** 
 * <p>This package was created using MATLAB Java Package. The classes included in this 
 * package are wrappers around MATLAB functions which were used when compiling this 
 * package in MATLAB. These classes have public methods that provide access to the 
 * MATLAB functions used by MATLAB Java Package during compilation.</p>
 * <h3><b>IMPORTANT : </b>What you need to use this package successfully :</h3>
 * <h3>MATLAB Runtime</h3>
 * <p>
 * <ul>
 * <li>MATLAB Runtime is the collection of platform specific native libraries required to execute M functions exported by the classes of this package</li>
 * <li>It can be made available either by installing MATLAB, MATLAB Compiler and MATLAB Java Package or by just running MCRInstaller executable</li>
 * <li>This package is compatible with MATLAB Runtime version 9.0 only.</li>
 * <li>Please contact the creator of this package for specific details about the MATLAB Runtime (e.g what version of MATLAB it originated with since MATLAB Runtime version is tied to the version of MATLAB)</li>
 * </ul> 
 * </p>
 * <p><b>NOTE: </b>Creating the first instance of one of the classes from this package is more time
 * consuming than creating subsequent instances since the native libraries from the MATLAB Runtime
 * must be loaded.</p>
 * <h3>javabuilder.jar</h3> 
 * <p>
 * <ul>
 * <li>Provides classes that act as the bridge between your application and the MATLAB Runtime</li>
 * <li>Located in the $MCR/toolbox/javabuilder/jar directory (where $MCR is the root of an installation of either MATLAB or MCRInstaller)</li>
 * <li>The <code>MSR</code> package will only work with the javabuilder.jar file included with MATLAB Runtime version 9.0</li>
 * </ul>
 * </p> 
 * <p><b>NOTE: </b><code>com.mathworks.toolbox.javabuilder.MWArray</code> is one of many data 
 * conversion classes provided in javabuilder.jar. MWArray is an abstract class representing a 
 * MATLAB array. Each MATLAB array type has a corresponding concrete class type in the 
 * MWArray class hierarchy. The public methods that represent MATLAB functions, for the classes 
 * in the <code>MSR</code> package, can take instances of these concrete classes as 
 * input. These methods can also take native Java primitive or array types as input. These native 
 * types are converted to appropriate MWArray types using data conversion rules provided by 
 * MATLAB Java Package. For instance, a Java primitive double is converted into an instance
 * of MWNumericArray (a subclass of MWArray).</p> 
 */
package MSR;
