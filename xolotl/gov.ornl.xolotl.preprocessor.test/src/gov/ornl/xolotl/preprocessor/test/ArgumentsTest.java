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
		final Arguments args;

		try {
			// Parse the empty string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {});
			
			// Check that the default maximum Helium cluster size is 8
			assertEquals(8, args.getMaxHeSize());
			
			// Check that the default maximum vacancy cluster size is 29
			assertEquals(29, args.getMaxVSize());
			
			// Check that the default maximum interstitial cluster size is 6
			assertEquals(6, args.getMaxISize());

			// Check the default material argument
			assertEquals("W100", args.getMaterial());

			// Check the default number of dimensions
			assertEquals("1", args.getDimensions());

			// Check if there is a startTemp argument
			assertEquals("1000", args.getStartTemp());

			// Check if there is a tempFile argument
			assertEquals(false, args.isTempFile());

			// Check if there is an heFlux argument
			assertEquals(false, args.isHeFlux());

			// Check if there is an fluxFile argument
			assertEquals(false, args.isFluxFile());

			// Check that the default perfHandler is std
			assertEquals("std", args.getPerfHandler());

			// Check if there is a vizHandler argument
			assertEquals("dummy", args.getVizHandler());

			// Check if there is a checkpoint argument
			assertEquals(false, args.isCheckpoint());

			// Check that the default networkFile is networkInit.h5
			assertEquals("networkInit.h5", args.getNetworkFile());

			// Check that the default stepSize is 1.0
			assertEquals("1.0", args.getStepSize());

			// Check the default petscArgs
			assertEquals(
					"-da_grid_x 10 -ts_final_time 50 -ts_dt 1.0e-12 "
							+ "-ts_max_steps 100 -ts_adapt_dt_max 10 -ts_max_snes_failures 200 "
							+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant "
							+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
					args.getPetscArgs());
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}
	}

	/**
	 * This operation tests that default parameter values are only overridden if
	 * they are specified via the command line and that the optional arguments 
	 * are only set if they are also specified
	 */
	@Test
	public void testSpecifiedArguments() {

		// Local Declarations
		final Arguments args;

		try {
			// Parse the specified string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {
					"--startTemp", "900", "--material", "W111", "--perfHandler",
					"dummy", "--dimensions", "2", "--maxHeSize", "7", "--maxVSize", 
					"30", "--maxISize", "5", "--checkpoint", "xolotlStop.h5", 
					"--stepSize", "3.0", "--initialV", "0.05"});
			
			// Check that the maximum Helium cluster size is 7
			assertEquals(7, args.getMaxHeSize());
			
			// Check that the maximum vacancy cluster size is 30
			assertEquals(30, args.getMaxVSize());
			
			// Check that the maximum interstitial cluster size is 5
			assertEquals(5, args.getMaxISize());
			
			// Check that the material is W111
			assertEquals("W111", args.getMaterial());
			
			// Check that the number of dimensions is 2
			assertEquals("2", args.getDimensions());
			
			// Check that the startTemp is 900
			assertEquals("900", args.getStartTemp());

			// Check if there is a tempFile argument
			assertEquals(false, args.isTempFile());

			// Check if there is an heFlux argument
			assertEquals(false, args.isHeFlux());

			// Check that the perfHandler is dummy
			assertEquals("dummy", args.getPerfHandler());

			// Check if there is a vizHandler argument
			assertEquals("dummy", args.getVizHandler());

			// Check if there is a checkpoint argument
			assertEquals(true, args.isCheckpoint());

			// Check the name of the file for the checkpoint
			assertEquals("xolotlStop.h5", args.getCheckpoint());

			// Check if there is an initial vacancy concentration argument
			assertEquals(true, args.isInitialV());

			// Check its value
			assertEquals("0.05", args.getInitialV());

			// Check that the default networkFile is networkInit.h5
			assertEquals("networkInit.h5", args.getNetworkFile());

			// Check that the default networkFile is networkInit.h5
			assertEquals("3.0", args.getStepSize());

			// Check the default petscArgs
			assertEquals(
					"-da_grid_x 10 -ts_final_time 50 -ts_dt 1.0e-12 "
							+ "-ts_max_steps 100 -ts_adapt_dt_max 10 -ts_max_snes_failures 200 "
							+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant "
							+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
					args.getPetscArgs());
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}
	}

}