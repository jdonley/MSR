/*
* MATLAB Compiler: 6.1 (R2015b)
* Date: Wed Oct 28 10:33:42 2015
* Arguments: "-B" "macro_default" "-W" "dotnet:MSR,Class1,0.0,private" "-T" "link:lib"
* "-d" "M:\MSR\MSR\for_testing" "-v" "class{Class1:M:\MSR\MSR.m}" "-a"
* "M:\MSR\+Orthogonal_Basis_Expansion\Gram_Schmidt.m" "-a"
* "M:\MSR\+Speaker_Setup\loudspeaker_setup.m" "-a"
* "M:\MSR\+Orthogonal_Basis_Expansion\multizone_soundfield_OBE.m" "-a"
* "M:\MSR\+Orthogonal_Basis_Expansion\spatial_zone.m" 
*/
using System;
using System.Reflection;
using System.IO;
using MathWorks.MATLAB.NET.Arrays;
using MathWorks.MATLAB.NET.Utility;

#if SHARED
[assembly: System.Reflection.AssemblyKeyFile(@"")]
#endif

namespace MSR
{

  /// <summary>
  /// The Class1 class provides a CLS compliant, MWArray interface to the MATLAB
  /// functions contained in the files:
  /// <newpara></newpara>
  /// M:\MSR\MSR.m
  /// </summary>
  /// <remarks>
  /// @Version 0.0
  /// </remarks>
  public class Class1 : IDisposable
  {
    #region Constructors

    /// <summary internal= "true">
    /// The static constructor instantiates and initializes the MATLAB Runtime instance.
    /// </summary>
    static Class1()
    {
      if (MWMCR.MCRAppInitialized)
      {
        try
        {
          Assembly assembly= Assembly.GetExecutingAssembly();

          string ctfFilePath= assembly.Location;

          int lastDelimiter= ctfFilePath.LastIndexOf(@"\");

          ctfFilePath= ctfFilePath.Remove(lastDelimiter, (ctfFilePath.Length - lastDelimiter));

          string ctfFileName = "MSR.ctf";

          Stream embeddedCtfStream = null;

          String[] resourceStrings = assembly.GetManifestResourceNames();

          foreach (String name in resourceStrings)
          {
            if (name.Contains(ctfFileName))
            {
              embeddedCtfStream = assembly.GetManifestResourceStream(name);
              break;
            }
          }
          mcr= new MWMCR("",
                         ctfFilePath, embeddedCtfStream, true);
        }
        catch(Exception ex)
        {
          ex_ = new Exception("MWArray assembly failed to be initialized", ex);
        }
      }
      else
      {
        ex_ = new ApplicationException("MWArray assembly could not be initialized");
      }
    }


    /// <summary>
    /// Constructs a new instance of the Class1 class.
    /// </summary>
    public Class1()
    {
      if(ex_ != null)
      {
        throw ex_;
      }
    }


    #endregion Constructors

    #region Finalize

    /// <summary internal= "true">
    /// Class destructor called by the CLR garbage collector.
    /// </summary>
    ~Class1()
    {
      Dispose(false);
    }


    /// <summary>
    /// Frees the native resources associated with this object
    /// </summary>
    public void Dispose()
    {
      Dispose(true);

      GC.SuppressFinalize(this);
    }


    /// <summary internal= "true">
    /// Internal dispose function
    /// </summary>
    protected virtual void Dispose(bool disposing)
    {
      if (!disposed)
      {
        disposed= true;

        if (disposing)
        {
          // Free managed resources;
        }

        // Free native resources
      }
    }


    #endregion Finalize

    #region Methods

    /// <summary>
    /// Provides a single output, 0-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR()
    {
      return mcr.EvaluateFunction("MSR", new MWArray[]{});
    }


    /// <summary>
    /// Provides a single output, 1-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res)
    {
      return mcr.EvaluateFunction("MSR", res);
    }


    /// <summary>
    /// Provides a single output, 2-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f)
    {
      return mcr.EvaluateFunction("MSR", res, f);
    }


    /// <summary>
    /// Provides a single output, 3-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R)
    {
      return mcr.EvaluateFunction("MSR", res, f, R);
    }


    /// <summary>
    /// Provides a single output, 4-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N);
    }


    /// <summary>
    /// Provides a single output, 5-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb);
    }


    /// <summary>
    /// Provides a single output, 6-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq);
    }


    /// <summary>
    /// Provides a single output, 7-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu);
    }


    /// <summary>
    /// Provides a single output, 8-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb);
    }


    /// <summary>
    /// Provides a single output, 9-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb);
    }


    /// <summary>
    /// Provides a single output, 10-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb);
    }


    /// <summary>
    /// Provides a single output, 11-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb);
    }


    /// <summary>
    /// Provides a single output, 12-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq);
    }


    /// <summary>
    /// Provides a single output, 13-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq);
    }


    /// <summary>
    /// Provides a single output, 14-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq);
    }


    /// <summary>
    /// Provides a single output, 15-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq, MWArray L)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L);
    }


    /// <summary>
    /// Provides a single output, 16-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq, MWArray L, MWArray Rl)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl);
    }


    /// <summary>
    /// Provides a single output, 17-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq, MWArray L, MWArray Rl, MWArray 
                 fmax)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax);
    }


    /// <summary>
    /// Provides a single output, 18-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <param name="phi">Input argument #18</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq, MWArray L, MWArray Rl, MWArray 
                 fmax, MWArray phi)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi);
    }


    /// <summary>
    /// Provides a single output, 19-input MWArrayinterface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <param name="phi">Input argument #18</param>
    /// <param name="phiL">Input argument #19</param>
    /// <returns>An MWArray containing the first output argument.</returns>
    ///
    public MWArray MSR(MWArray res, MWArray f, MWArray R, MWArray N, MWArray Wb, MWArray 
                 Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray disb, MWArray PWangb, 
                 MWArray rq, MWArray angq, MWArray disq, MWArray L, MWArray Rl, MWArray 
                 fmax, MWArray phi, MWArray phiL)
    {
      return mcr.EvaluateFunction("MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi, phiL);
    }


    /// <summary>
    /// Provides the standard 0-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", new MWArray[]{});
    }


    /// <summary>
    /// Provides the standard 1-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res);
    }


    /// <summary>
    /// Provides the standard 2-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f);
    }


    /// <summary>
    /// Provides the standard 3-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R);
    }


    /// <summary>
    /// Provides the standard 4-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N);
    }


    /// <summary>
    /// Provides the standard 5-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb);
    }


    /// <summary>
    /// Provides the standard 6-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq);
    }


    /// <summary>
    /// Provides the standard 7-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu);
    }


    /// <summary>
    /// Provides the standard 8-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb);
    }


    /// <summary>
    /// Provides the standard 9-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb);
    }


    /// <summary>
    /// Provides the standard 10-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb);
    }


    /// <summary>
    /// Provides the standard 11-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb);
    }


    /// <summary>
    /// Provides the standard 12-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq);
    }


    /// <summary>
    /// Provides the standard 13-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq);
    }


    /// <summary>
    /// Provides the standard 14-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq);
    }


    /// <summary>
    /// Provides the standard 15-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq, MWArray 
                   L)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L);
    }


    /// <summary>
    /// Provides the standard 16-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq, MWArray 
                   L, MWArray Rl)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl);
    }


    /// <summary>
    /// Provides the standard 17-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq, MWArray 
                   L, MWArray Rl, MWArray fmax)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax);
    }


    /// <summary>
    /// Provides the standard 18-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <param name="phi">Input argument #18</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq, MWArray 
                   L, MWArray Rl, MWArray fmax, MWArray phi)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi);
    }


    /// <summary>
    /// Provides the standard 19-input MWArray interface to the MSR MATLAB function.
    /// </summary>
    /// <remarks>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return.</param>
    /// <param name="res">Input argument #1</param>
    /// <param name="f">Input argument #2</param>
    /// <param name="R">Input argument #3</param>
    /// <param name="N">Input argument #4</param>
    /// <param name="Wb">Input argument #5</param>
    /// <param name="Wq">Input argument #6</param>
    /// <param name="Wu">Input argument #7</param>
    /// <param name="rb">Input argument #8</param>
    /// <param name="angb">Input argument #9</param>
    /// <param name="disb">Input argument #10</param>
    /// <param name="PWangb">Input argument #11</param>
    /// <param name="rq">Input argument #12</param>
    /// <param name="angq">Input argument #13</param>
    /// <param name="disq">Input argument #14</param>
    /// <param name="L">Input argument #15</param>
    /// <param name="Rl">Input argument #16</param>
    /// <param name="fmax">Input argument #17</param>
    /// <param name="phi">Input argument #18</param>
    /// <param name="phiL">Input argument #19</param>
    /// <returns>An Array of length "numArgsOut" containing the output
    /// arguments.</returns>
    ///
    public MWArray[] MSR(int numArgsOut, MWArray res, MWArray f, MWArray R, MWArray N, 
                   MWArray Wb, MWArray Wq, MWArray Wu, MWArray rb, MWArray angb, MWArray 
                   disb, MWArray PWangb, MWArray rq, MWArray angq, MWArray disq, MWArray 
                   L, MWArray Rl, MWArray fmax, MWArray phi, MWArray phiL)
    {
      return mcr.EvaluateFunction(numArgsOut, "MSR", res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi, phiL);
    }


    /// <summary>
    /// Provides an interface for the MSR function in which the input and output
    /// arguments are specified as an array of MWArrays.
    /// </summary>
    /// <remarks>
    /// This method will allocate and return by reference the output argument
    /// array.<newpara></newpara>
    /// M-Documentation:
    /// </remarks>
    /// <param name="numArgsOut">The number of output arguments to return</param>
    /// <param name= "argsOut">Array of MWArray output arguments</param>
    /// <param name= "argsIn">Array of MWArray input arguments</param>
    ///
    public void MSR(int numArgsOut, ref MWArray[] argsOut, MWArray[] argsIn)
    {
      mcr.EvaluateFunction("MSR", numArgsOut, ref argsOut, argsIn);
    }



    /// <summary>
    /// This method will cause a MATLAB figure window to behave as a modal dialog box.
    /// The method will not return until all the figure windows associated with this
    /// component have been closed.
    /// </summary>
    /// <remarks>
    /// An application should only call this method when required to keep the
    /// MATLAB figure window from disappearing.  Other techniques, such as calling
    /// Console.ReadLine() from the application should be considered where
    /// possible.</remarks>
    ///
    public void WaitForFiguresToDie()
    {
      mcr.WaitForFiguresToDie();
    }



    #endregion Methods

    #region Class Members

    private static MWMCR mcr= null;

    private static Exception ex_= null;

    private bool disposed= false;

    #endregion Class Members
  }
}
