# Copyright (c) 2009-2016 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

""" HOOMD-blue python API

:py:mod:`hoomd` provides a high level user interface for executing
simulations using HOOMD::

    import hoomd
    from hoomd import md
    hoomd.context.initialize("")
    hoomd.init.create_random(N=100, phi_p=0.1)
    lj = md.pair.lj(r_cut=2.5)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
    hoomd.md.integrate.mode_standard(dt=0.005)
    hoomd.md.integrate.nvt(group=hoomd.group.all(), T=1.2, tau=0.5)
    hoomd.run(100)

.. rubric:: Stability

:py:mod:`hoomd` is **stable**. When upgrading from version 2.x to 2.y (y > x),
existing job scripts that follow *documented* interfaces for functions and classes
will not require any modifications. **Maintainer:** Joshua A. Anderson

.. attention::

    This stability guaruntee only applies to modules in the :py:mod:`hoomd` package.
    Subpackages (:py:mod:`hoomd.hpmc`, :py:mod:`hoomd.md`, etc...) may or may not
    have a stable API. The documentation for each subpackage specifies the level of
    API stability it provides.
"""

# Maintainer: joaander
import sys;
import ctypes;
import os;

# need to import HOOMD with RTLD_GLOBAL in python sitedir builds
if not ('NOT_HOOMD_PYTHON_SITEDIR' in os.environ):
    flags = sys.getdlopenflags();
    sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL);

from hoomd import _hoomd;

if not ('NOT_HOOMD_PYTHON_SITEDIR' in os.environ):
    sys.setdlopenflags(flags);

from hoomd import meta
from hoomd import context
from hoomd import cite
from hoomd import analyze
from hoomd import benchmark
from hoomd import comm
from hoomd import compute
from hoomd import data
from hoomd import dump
from hoomd import group
from hoomd import init
from hoomd import integrate
from hoomd import option
from hoomd import update
from hoomd import util
from hoomd import variant
from hoomd import lattice

# from hoomd._hoomd import WalltimeLimitReached; TODO: adios_boost, re-enable this with custom exception code

# output the version info on import
context.msg.notice(1, _hoomd.output_version_info())

# ensure creation of global bibliography to print HOOMD base citations
cite._ensure_global_bib()

_default_excepthook = sys.excepthook;

## \internal
# \brief Override pythons except hook to abort MPI runs
def _hoomd_sys_excepthook(type, value, traceback):
    _default_excepthook(type, value, traceback);
    sys.stderr.flush();
    if context.exec_conf is not None:
        _hoomd.abort_mpi(context.exec_conf);

__version__ = "{0}.{1}.{2}".format(*_hoomd.__version__)

def run(tsteps, profile=False, limit_hours=None, limit_multiple=1, callback_period=0, callback=None, quiet=False):
    """ Runs the simulation for a given number of time steps.

    Args:

        tsteps (int): Number of time steps to advance the simulation.
        profile (bool): Set to True to enable high level profiling output at the end of the run.
        profile limit_hours (float): If not None, limit this run to a given number of hours.
        limit_multiple (int): When stopping the run due to walltime limits, only stop when the time step is a
                              multiple of limit_multiple.
        callback (python callable): Sets a Python function to be called regularly during a run.
        callback_period (int): Sets the period, in time steps, between calls made to ``callback``.
        quiet (bool): Set to True to disable the status information printed to the screen by the run.

    Example::

            hoomd.run(10)
            hoomd.run(10e6, limit_hours=1.0/3600.0, limit_multiple=10)
            hoomd.run(10, profile=True)
            hoomd.run(10, quiet=True)
            hoomd.run(10, callback_period=2, callback=lambda step: print(step))

    Execute the :py:func:`run()` command to advance the simulation forward in time.
    During the run, all previously specified analyzers, updaters and the integrator
    are executed at the specified regular periods.

    After :py:func:`run()` completes, you may change parameters of the simulation
    and continue the simulation by executing :py:func:`run()` again. Time steps are added
    cumulatively, so calling ``run(1000)`` and then ``run(2000)`` would run the simulation
    up to time step 3000.

    :py:func:`run()` cannot be executed before the system is initialized. In most
    cases, :py:func:`run()` should only be called after after pair forces, bond forces,
    and an integrator are specified.

    When ``profile`` is **True**, a detailed breakdown of how much time was spent in each
    portion of the calculation is printed at the end of the run. Collecting this timing information
    slows the simulation.

    **Wallclock limited runs:**

    There are a number of mechanisms to limit the time of a running hoomd script. Use these in a job
    queuing environment to allow your script to cleanly exit before reaching the system enforced walltime limit.

    Force :py:func:`run()` to end only on time steps that are a multiple of ``limit_mulitple``. Set this to the period at which you
    dump restart files so that you always end a :py:func:`run()` cleanly at a point where you can restart from. Use
    ``phase=0`` on logs, file dumps, and other periodic tasks. With ``phase=0``, these tasks will continue on the same
    sequence regardless of the restart period.

    Set the environment variable ``HOOMD_WALLTIME_STOP`` prior to starting a hoomd script to stop the :py:func:`run()` at a given wall
    clock time. :py:func:`run()` monitors performance and tries to ensure that it will end *before* ``HOOMD_WALLTIME_STOP``. This
    environment variable works even with multiple stages of runs in a script (use :py:func:`run_upto()`. Set the variable to
    a unix epoch time. For example in a job script that should run 12 hours, set ``HOOMD_WALLTIME_STOP`` to 12 hours from
    now, minus 10 minutes to allow for job cleanup::

        export HOOMD_WALLTIME_STOP=$((`date +%s` + 12 * 3600 - 10 * 60))

    When using ``HOOMD_WALLTIME_STOP``, :py:func:`run()` will throw the exception ``WalltimeLimitReached`` if it exits due to the walltime
    limit.

    ``limit_hours`` is another way to limit the length of a :py:func:`run()`. Set it to a number of hours (use fractional values for
    minutes) to limit this particular :py:func:`run()` to that length of time. This is less useful than ``HOOMD_WALLTIME_STOP`` in a
    job queuing environment.

    **Callbacks:**

    If ``callback`` is set to a Python function then this function will be called regularly
    at ``callback_period`` intervals. The callback function must receive one integer as argument
    and can return an integer. The argument passed to the callback is the current time step number.
    If the callback function returns a negative number, the run is immediately aborted.

    If ``callback_period`` is set to 0 (the default) then the callback is only called
    once at the end of the run. Otherwise the callback is executed whenever the current
    time step number is a multiple of ``callback_period``.
    """

    if not quiet:
        util.print_status_line();
    # check if initialization has occured
    if not init.is_initialized():
        context.msg.error("Cannot run before initialization\n");
        raise RuntimeError('Error running');

    if context.current.integrator is None:
        context.msg.warning("Starting a run without an integrator set");
    else:
        context.current.integrator.update_forces();
        context.current.integrator.update_methods();
        context.current.integrator.update_thermos();

    # update autotuner parameters
    context.current.system.setAutotunerParams(context.options.autotuner_enable, int(context.options.autotuner_period));

    for logger in context.current.loggers:
        logger.update_quantities();
    context.current.system.enableProfiler(profile);
    context.current.system.enableQuietRun(quiet);

    # update all user-defined neighbor lists
    for nl in context.current.neighbor_lists:
        nl.update_rcut()
        nl.update_exclusions_defaults()

    # detect 0 hours remaining properly
    if limit_hours == 0.0:
        context.msg.warning("Requesting a run() with a 0 time limit, doing nothing.\n");
        return;
    if limit_hours is None:
        limit_hours = 0.0

    if not quiet:
        context.msg.notice(1, "** starting run **\n");
    context.current.system.run(int(tsteps), callback_period, callback, limit_hours, int(limit_multiple));
    if not quiet:
        context.msg.notice(1, "** run complete **\n");

def run_upto(step, **keywords):
    """Runs the simulation up to a given time step number.


    Args:

        step (int): Final time step of the simulation which to run
        keywords (see below): Catch for all keyword arguments to pass on to :py:func:`run()`

    :py:func:`run_upto()` runs the simulation, but only until it reaches the given time step. If the simulation has already
    reached the specified step, a message is printed and no simulation steps are run.

    It accepts all keyword options that :py:func:`run()` does.

    Examples::

        run_upto(1000)
        run_upto(10000, profile=True)
        run_upto(1e9, limit_hours=11)
    """
    if 'quiet' in keywords and not keywords['quiet']:
        util.print_status_line();
    # check if initialization has occured
    if not init.is_initialized():
        context.msg.error("Cannot run before initialization\n");
        raise RuntimeError('Error running');

    # determine the number of steps to run
    step = int(step);
    cur_step = context.current.system.getCurrentTimeStep();

    if cur_step >= step:
        context.msg.notice(2, "Requesting run up to a time step that has already passed, doing nothing\n");
        return;

    n_steps = step - cur_step;

    util.quiet_status();
    run(n_steps, **keywords);
    util.unquiet_status();

def get_step():
    """ Get the current simulation time step.

    Returns:
        The current simulation time step.

    Example::

            print(hoomd.get_step())
    """

    # check if initialization has occurred
    if not init.is_initialized():
        context.msg.error("Cannot get step before initialization\n");
        raise RuntimeError('Error getting step');

    return context.current.system.getCurrentTimeStep();

# Check to see if we are built without MPI support and the user used mpirun
if (not _hoomd.is_MPI_available()) and ('OMPI_COMM_WORLD_RANK' in os.environ or 'MV2_COMM_WORLD_LOCAL_RANK' in os.environ):
    print('HOOMD-blue is built without MPI support, but seems to have been launched with mpirun');
    print('exiting now to prevent many sequential jobs from starting');
    raise RuntimeError('Error launching hoomd')
