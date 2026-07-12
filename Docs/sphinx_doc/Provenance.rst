.. _sec:Provenance:

******************************
Provenance and Restart Lineage
******************************

Purpose
=======

ERF records the identity of each native checkpoint and primary native AMReX
plotfile. The record also describes the ERF invocation that wrote the output
and the known restart lineage. Provenance is auxiliary metadata. It does not
change model state or restart calculations.

Identifiers
===========

ERF uses three identities:

* A **simulation UUID** identifies one restart lineage.
* An **execution UUID** identifies one ERF invocation.
* An **artifact UUID** identifies one checkpoint or plotfile directory.

The parent execution UUID identifies the immediate producer of a restart.
The complete known execution lineage records every known execution in order.
Both fields are present because one execution can write several checkpoints.
The source checkpoint UUID identifies the exact checkpoint used for restart.

The following fields appear in every provenance block.

.. list-table:: Provenance fields
   :header-rows: 1
   :widths: 25 50 25

   * - Field
     - Meaning
     - Lifetime
   * - simulation UUID
     - Root identity shared by one restart lineage
     - Entire lineage
   * - execution UUID
     - Identity of one ERF invocation
     - One invocation
   * - parent execution UUID
     - Immediate producer execution
     - One restart edge
   * - execution lineage
     - Ordered known execution history
     - Accumulates
   * - artifact UUID
     - Identity of one checkpoint or plotfile
     - One artifact
   * - source checkpoint UUID
     - Exact checkpoint used for restart
     - One invocation
   * - restart generation
     - Number of restarts from the root, or ``-1`` if unknown
     - One invocation
   * - lineage status
     - Completeness or the reason metadata could not be recovered
     - One invocation

Cold Starts
===========

On a cold start, the simulation UUID and execution UUID are equal. The
lineage contains one UUID, the parent is empty, and the restart generation is
zero. The lineage status is ``complete``.

Checkpoint Restarts
===================

Suppose execution ``A`` writes checkpoint ``CHK_A`` and execution ``B``
restarts from it. Execution ``B`` keeps the simulation UUID from ``A``, sets
its parent execution UUID to ``A``, appends ``B`` to the lineage, and stores
``CHK_A`` as its source checkpoint UUID. Its restart generation is one.

Repeated Restarts
=================

The full known lineage is copied into every new artifact. A sequence
``A -> B -> C`` therefore appears as ``[A, B, C]`` in artifacts written by
``C``. ERF does not need the ``A`` directory to inspect ``C``'s lineage.

Branching Restarts
==================

Two executions may restart from the same checkpoint. For example, children
``C`` and ``D`` of ``B`` share the simulation UUID, parent ``B``, and known
lineage prefix ``[A, B]``. They have different execution UUIDs and their
lineages end in ``C`` and ``D`` respectively.

Legacy Checkpoints
==================

Older checkpoints have no provenance block. Missing ``job_info``, a missing
block, malformed fields, invalid UUIDs, failed lineage validation, an
unsupported schema, or a non-checkpoint artifact record produces one warning
on the I/O rank. ERF then reads the physical checkpoint normally.

The new execution starts a known lineage segment. Its lineage contains only
the current execution, its restart generation is ``-1``, and its lineage is
incomplete. The status records the specific failure. A later restart from an
artifact written by this execution preserves the known segment and appends
the child, while keeping the generation unknown.

Native File Locations
=====================

The block is embedded near the top of the existing human-readable
``job_info`` file. ERF writes it to native checkpoints and to primary native
2D and 3D AMReX plotfiles. It does not write a sidecar.

NetCDF output, NetCDF column files, 2D diagnostic JSON, subvolumes, staggered
velocity plotfiles, boundary-plane files, and particle-only files are outside
this implementation.

Schema and Validation
=====================

Schema version one uses fixed-order ``key=value`` lines between
``ERF_PROVENANCE_BEGIN`` and ``ERF_PROVENANCE_END``. Paths are split at the
first equals sign. The parser accepts CRLF line endings, ignores unknown keys,
rejects duplicate required keys, and validates UUIDs and lineage invariants.
The checkpoint ``Header`` remains unchanged.

UUIDs are generated on the I/O rank with an independent host-side generator
and broadcast to all ranks. This generator does not use AMReX's scientific
random-number stream. UUIDs identify records but do not prove authenticity;
they are not signatures or cryptographic evidence.

Scope and Limits
================

This feature deliberately does not compute configuration hashes, classify
configuration changes, digest input files or executables, export W3C PROV or
RO-Crate records, sign metadata, or reject a restart based on provenance.
Configuration identity and restart lineage remain separate concerns.
