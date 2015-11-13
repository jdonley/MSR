/*
 * MATLAB Compiler: 6.1 (R2015b)
 * Date: Tue Oct 27 18:52:12 2015
 * Arguments: "-B" "macro_default" "-W" "java:MSR,Class1" "-T" "link:lib" "-d" 
 * "M:\\MSR\\MSR\\for_testing" "class{Class1:M:\\MSR\\MSR.m}" "-a" 
 * "M:\\MSR\\+Orthogonal_Basis_Expansion\\Gram_Schmidt.m" "-a" 
 * "M:\\MSR\\+Speaker_Setup\\loudspeaker_setup.m" "-a" 
 * "M:\\MSR\\+Orthogonal_Basis_Expansion\\multizone_soundfield_OBE.m" "-a" 
 * "M:\\MSR\\+Orthogonal_Basis_Expansion\\spatial_zone.m" 
 */

package MSR;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class MSRMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "MSR_51226730D325D8B0D8D9D0723E8426D4";
    
    /** Component name */
    private static final String sComponentName = "MSR";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(MSRMCRFactory.class)
        );
    
    
    private MSRMCRFactory()
    {
        // Never called.
    }
    
    public static MWMCR newInstance(MWComponentOptions componentOptions) throws MWException
    {
        if (null == componentOptions.getCtfSource()) {
            componentOptions = new MWComponentOptions(componentOptions);
            componentOptions.setCtfSource(sDefaultComponentOptions.getCtfSource());
        }
        return MWMCR.newInstance(
            componentOptions, 
            MSRMCRFactory.class, 
            sComponentName, 
            sComponentId,
            new int[]{9,0,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
