package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import org.junit.Test;

import gov.ornl.xolotl.preprocessor.Arguments;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class is responsible for testing the Arguments command line interface.
 */
public class ArgumentsTest {

	/**
	 * This operation checks that the default argument values are used if no
	 * command line options are specified.
	 */
	@Test
	public void testDefaultArguments() {
		// Local Declarations
		Arguments args;

		try {
			// Parse the empty string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {});
			
			// Check that the default maximum Helium cluster size is 8
			assertEquals(8, args.getMaxHeSize());
			
			// Check that the default maximum vacancy cluster size is 29
			assertEquals(29, args.getMaxVSize());
			
			// Check that the default maximum interstitial cluster size is 6
			assertEquals(6, args.getMaxISize());

			// Check if there is a phase cut argument
			assertEquals(false, args.isPhaseCut());

			// Check if there is a startTemp argument
			assertEquals("1000", args.getStartTemp());

			// Check that the default perfHandler is std
			assertEquals("std", args.getPerfHandler());

			// Check if there is a vizHandler argument
			assertEquals("dummy", args.getVizHandler());

			// Check the default petscArgs
			assertEquals("-ts_final_time 1.0 -ts_dt 1.0e-12 "
					+ "-ts_max_steps 100 -ts_adapt_dt_max 1.0e-6 -ts_max_snes_failures 200 "
					+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type sor "
					+ "-fieldsplit_1_pc_type redundant -ts_monitor",
					args.getPetscArgs());

			// Check that the default networkFile is networkInit.h5
			assertEquals("networkInit.h5", args.getNetworkFile());

			// Check the default number of dimensions
			assertEquals("1", args.getDimensions());

			// Check the default number of grid points in the x direction
			assertEquals(20, args.getNxGrid());

			// Check the default number of grid points in the y direction
			assertEquals(0, args.getNyGrid());

			// Check the default number of grid points in the z direction
			assertEquals(0, args.getNzGrid());

			// Check that the default xStepSize is 1.0
			assertEquals(1.0, args.getXStepSize(), 0.001);

			// Check that the default yStepSize is 0.0
			assertEquals(0.0, args.getYStepSize(), 0.001);

			// Check that the default zStepSize is 0.0
			assertEquals(0.0, args.getZStepSize(), 0.001);

			// Check the default material argument
			assertEquals("W100", args.getMaterial());

			// Check the default flux argument
			assertEquals("4.0e7", args.getFlux());

			// Check if there is a tempFile argument
			assertEquals(false, args.isTempFile());

			// Check if there is an fluxFile argument
			assertEquals(false, args.isFluxFile());

			// Check if there is a checkpoint argument
			assertEquals(false, args.isCheckpoint());

			// Check if there is a initial vacancy concentration argument
			assertEquals(false, args.isInitialV());
						
			// Check the default Td argument
			assertEquals("80", args.getThresholdEnergy());
			
			// Check the default krypton fluence argument
			assertEquals("0", args.getKrFluence());
		
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		} 
		
		return;
	}

	/**
	 * This operation tests that default parameter values are only overridden if
	 * they are specified via the command line and that the optional arguments 
	 * are only set if they are also specified
	 */
	@Test
	public void testSpecifiedArguments() {
		// Local Declarations
		Arguments args;

		try {
			// Parse the specified string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {
				"--maxHeSize", "7", "--maxVSize", "30", "--maxISize", "5", "--phaseCut", 
				"--startTemp", "900", "--perfHandler", "dummy", "--vizHandler", "std", 
				"--petscArgs=-plot", "--networkFile", "net.h5",
				"--dimensions", "2", "--nxGrid", "50", "--nyGrid", "10", "--nzGrid", "30", 
				"--xStepSize", "0.2", "--yStepSize", "1.5", "--zStepSize", "10.0", 
				"--material", "W111", "--tempFile", "temp.dat", "--flux", "5.0e5", 
				"--fluxFile", "flux.dat", "--checkpoint", "xolotlStop.h5", 
				"--initialV", "0.05", "--krFluence", "8.0e15", "--thresholdEnergy", "120" });
			
			// Check that the maximum Helium cluster size is 7
			assertEquals(7, args.getMaxHeSize());
			
			// Check that the maximum vacancy cluster size is 30
			assertEquals(30, args.getMaxVSize());
			
			// Check that the maximum interstitial cluster size is 5
			assertEquals(5, args.getMaxISize());

			// Check that the phase cut method is activated
			assertEquals(true, args.isPhaseCut());

			// Check that the temperature is 900
			assertEquals("900", args.getStartTemp());

			// Check that the perfHandler is dummy
			assertEquals("dummy", args.getPerfHandler());

			// Check that the vizHandler is std
			assertEquals("std", args.getVizHandler());

			// Check the petscArgs
			assertEquals("-plot", args.getPetscArgs());

			// Check that the networkFile is net.h5
			assertEquals("net.h5", args.getNetworkFile());
			
			// Check that the number of dimensions is 2
			assertEquals("2", args.getDimensions());
			
			// Check that the number of grid points in the x direction is 50
			assertEquals(50, args.getNxGrid());
			
			// Check that the number of grid points in the y direction is 10
			assertEquals(10, args.getNyGrid());
			
			// Check that the number of grid points in the z direction is 30
			assertEquals(30, args.getNzGrid());

			// Check that the xStepSize was set to 0.2
			assertEquals(0.2, args.getXStepSize(), 0.001);

			// Check that the yStepSize was set to 1.5
			assertEquals(1.5, args.getYStepSize(), 0.001);

			// Check that the zStepSize was set to 10.0
			assertEquals(10.0, args.getZStepSize(), 0.001);
		
			// Check that the material is W111
			assertEquals("W111", args.getMaterial());

			// Check that the flux argument is 5.0e5
			assertEquals("5.0e5", args.getFlux());

			// Check if there is a tempFile argument
			assertEquals(true, args.isTempFile());

			// Check that the tempFile argument is temp.dat
			assertEquals("temp.dat", args.getTempFile());

			// Check if there is an fluxFile argument
			assertEquals(true, args.isFluxFile());

			// Check that the fluxFile argument is flux.dat
			assertEquals("flux.dat", args.getFluxFile());
		
			// Check if there is a checkpoint argument
			assertEquals(true, args.isCheckpoint());

			// Check the name of the file for the checkpoint
			assertEquals("xolotlStop.h5", args.getCheckpoint());
			
			// Check if there is an initial vacancy concentration argument
			assertEquals(true, args.isInitialV());

			// Check its value
			assertEquals("0.05", args.getInitialV());

			// Check the krypton fluence is set to 8.0e15
			assertEquals("8.0e15", args.getKrFluence());

			// Check the threshold displacement energy is set to 120
			assertEquals("120", args.getThresholdEnergy());
				
		} 
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		} 
		
		return;
	}
}