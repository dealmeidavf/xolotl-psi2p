package gov.ornl.xolotl.preprocessor;

import uk.co.flamingpenguin.jewel.cli.Option;

/**
 * This interface creates command line options for the parameters that are needed 
 * to run Xolotl.  If the options are not specified via the command line, then 
 * the default values are used.
 * 
 * Each Option is constructed with a name (longName), a default value, and a help description.
 * E.g.  
 * 		@Option(longName="heFlux", defaultValue = "2.5e27",
 *				description = "The Helium insertion flux")
 * The above example creates the command line option '--heFlux'.  If this option is
 * not specified, then the default value is used.
 *
 */
public interface Arguments {
	
	@Option(longName="material", defaultValue = "W",
			description = "{W, Fe}  Which material will be used, Tungsten or Iron")
	String getMaterial();
	
	@Option(longName="startTemp", defaultValue = "1000",
			description = "The temperature (in Kelvin) will be the constant"
			+ "value specified")
	String getStartTemp();
	
	// The tempFile option is optional and should only be used when the user wishes 
	// to pass a file containing a temperature profile to Xolotl
	@Option(longName="tempFile", description = "<tempFileName>  A temperature profile is given "
			+ "by the specified file, then linear interpolation is used to fit the data")
	String getTempFile();
	boolean isTempFile();
	
	@Option(longName="heFlux", defaultValue = "2.5e27",
			description = "The Helium insertion flux")
	String getHeFlux();
	
	@Option(longName="heFluence", defaultValue = "1.0e19",
			description = "The Helium fluence value that will be used")
	String getHeFluence();
	
	@Option(longName="perfHandler", defaultValue = "dummy",
			description = "{dummy, std}  Which set of handlers to use")
	String getPerfHandler();
	
	@Option(longName="vizHandler", defaultValue = "dummy",
			description = "{dummy, std}  Which set of handlers to use")
	String getVizHandler();
	
	@Option(longName="petscArgs", defaultValue = "-da_grid_x 10 -ts_final_time 1000"
			+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200"
			+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant"
			+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
			description = "List of arguments to be passed to PETSc")
	String getPetscArgs();
	
	@Option(longName="networkFile", defaultValue = "test.h5",
			description = "<fileName>  The file containing the network")
	String getNetworkFile();
	
	@Option(longName="checkpoint", defaultValue = "true",
			description = "{true, false}  Write the output files")
	String getCheckpoint();
	boolean isCheckpoint();
	
	// Print usage help
	@Option(helpRequest = true, description = "display help", shortName = "h")
	boolean getHelp();
	
}
