# -*- coding: iso-8859-1 -*-

from hoomd import *
from hoomd import md
context.initialize()
import unittest
import os


r_colloid = 1

# tests for md.update.dynamic_bond
class update_dynamic_bond (unittest.TestCase):
    def setUp(self):
        snap = data.make_snapshot(N=2, box=data.boxdim(Lx=80, Ly=80, Lz=80),
                                  bond_types=['polymer'], particle_types=['colloid'])

        if comm.get_rank() == 0:
            snap.particles.diameter[:] = [r_colloid*2]*2
            snap.particles.position[0] = [-1, 0, 0]
            snap.particles.position[1] = [1, 0, 0]

        self.s = init.read_snapshot(snap);
        self.nl = md.nlist.tree();

    # tests basic creation of the updater
    def test_create(self):
        # nl = md.nlist.tree();
        dybond = md.update.dynamic_bond(group=group.all(), nlist=self.nl, seed=1, period=1)
        run(10);

    # tests formation of a bond within a cutoff radius
    def test_formation(self):
        dybond = md.update.dynamic_bond(group.all(), nlist=self.nl, seed=1, period=1)
        dybond.set_params(r_cut=2.0, bond_type='harmonic', prob_form=1, prob_break=0)

    # tests breakage of a bond within a cutoff radius
    def test_breakage(self):
        dybond = md.update.dynamic_bond(group.all(), nlist=self.nl, seed=1, period=1)
        dybond.set_params(r_cut=2.0, bond_type='harmonic', prob_form=0, prob_break=1)

    # test that no bonds are formed outside the rcut
    def test_outside_cutoff(self):
        dybond = md.update.dynamic_bond(group.all(), nlist=self.nl, seed=1, period=1)
        dybond.set_params(r_cut=1.0, bond_type='harmonic', prob_form=1, prob_break=0)

    def tearDown(self):
        del self.s
        context.initialize();

if __name__ == '__main__':
    unittest.main(argv = ['test.py', '-v'])



