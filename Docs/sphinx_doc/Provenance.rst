.. _sec:Provenance:

******************************
Provenance and Restart Lineage
******************************

Purpose
=======

ERF records the identity of each native checkpoint and primary native AMReX
plotfile. The record describes the ERF invocation that wrote the output and
the known ancestry path that led to it. Provenance is auxiliary metadata. It
does not change model state or restart calculations.

Identifiers
===========

The **simulation UUID** identifies one rooted restart family. A family may
contain several branches. The **execution UUID** identifies one ERF process
invocation. The **artifact UUID** identifies one checkpoint or primary native
plotfile directory.

The **parent execution UUID** identifies the execution that wrote the source
checkpoint used to start the current execution. The **execution lineage** is
the ordered known ancestry path that leads to the current execution. The
**source checkpoint UUID** identifies the exact checkpoint artifact used to
start it. The parent field gives the immediate graph edge; the lineage makes
the known ancestry path self-contained in every artifact. The source UUID is
separate because one execution can write several checkpoints.

Cold Starts
===========

On a cold start, the simulation UUID and execution UUID are equal. The
lineage contains one UUID, the parent is empty, the restart generation is
zero, and the lineage status is ``complete``.

Checkpoint Restarts
===================

Suppose execution ``A`` writes checkpoint ``CHK_A`` and execution ``B``
restarts from it. Execution ``B`` keeps the simulation UUID from ``A``, sets
its parent to ``A``, appends ``B`` to its ancestry path, and stores ``CHK_A``
as its source checkpoint UUID. Its restart generation is one.

Repeated Restarts
=================

The known ancestry path is copied into every new artifact. A sequence
``A -> B -> C`` therefore appears as ``[A, B, C]`` in artifacts written by
``C``. ERF does not need the ``A`` directory to inspect ``C``'s ancestry.

Branching Restarts
==================

A restart family can branch:

.. code-block:: text

   A -> B -> C
          \
           -> D

Both branches share the same simulation UUID. The ancestry path for ``C`` is
``[A, B, C]``. The path for ``D`` is ``[A, B, D]``.

Persisted Schema
================

The writer emits schema-1 keys in a stable order. The parser identifies
fields by key, so a valid block may use another order. Marker lines must
equal the complete marker after an optional trailing carriage return. The
parser splits values at the first equals sign, ignores unknown keys, rejects
duplicate required keys, and validates typed values and lineage invariants.

Every schema-1 key is listed below.

.. list-table:: Schema-1 provenance keys
   :header-rows: 1
   :widths: 30 45 25

   * - Key
     - Meaning
     - Scope or notes
   * - ``schema_version``
     - Version of the embedded provenance schema.
     - Integer; this implementation writes ``1``.
   * - ``simulation_uuid``
     - Identity of the rooted restart family.
     - Stable across complete restarts and branches.
   * - ``execution_uuid``
     - Identity of one ERF process invocation.
     - New for each invocation.
   * - ``parent_execution_uuid``
     - Execution that wrote the source checkpoint.
     - Empty on a cold start or unavailable parent metadata.
   * - ``execution_lineage``
     - Ordered known ancestry path ending at ``execution_uuid``.
     - Comma-separated UUIDs in the file.
   * - ``restart_generation``
     - Number of restart edges from the known root.
     - ``-1`` when the absolute generation is unknown.
   * - ``lineage_complete``
     - Whether the path begins at the original cold start.
     - ``true`` or ``false``.
   * - ``lineage_status``
     - Completeness state or recovery-failure reason.
     - Stable lowercase token.
   * - ``source_checkpoint_uuid``
     - UUID of the exact checkpoint used for restart.
     - Empty on a cold start or when unavailable.
   * - ``source_checkpoint_path``
     - Configured path of the restart checkpoint.
     - Location metadata; not content identity.
   * - ``execution_start_utc``
     - UTC start time of the ERF invocation.
     - ISO 8601 format.
   * - ``artifact_uuid``
     - Identity of this output artifact.
     - New for every checkpoint or primary native plotfile.
   * - ``artifact_type``
     - Type of this output artifact.
     - ``checkpoint``, ``plotfile_2d``, or ``plotfile_3d``.
   * - ``artifact_step``
     - Level-0 step associated with the output.
     - Nonnegative integer.
   * - ``artifact_time_seconds``
     - Simulation time associated with the output.
     - Serialized with stable double precision.
   * - ``artifact_created_utc``
     - UTC time when ERF finalized the artifact metadata.
     - ISO 8601 format.

Status Values
=============

The ``lineage_status`` value uses these stable tokens:

.. list-table:: Lineage status tokens
   :header-rows: 1
   :widths: 32 68

   * - Token
     - Meaning
   * - ``complete``
     - Ancestry is known back to the cold start.
   * - ``incomplete_ancestor``
     - The parent record was valid, but its earlier ancestry was incomplete.
   * - ``missing_job_info``
     - The source checkpoint had no readable ``job_info``.
   * - ``missing_provenance_block``
     - ``job_info`` existed but had no exact provenance block.
   * - ``malformed_provenance``
     - The block could not be parsed or failed validation.
   * - ``unsupported_schema``
     - The block used an unsupported schema version.
   * - ``artifact_type_mismatch``
     - The record did not describe a checkpoint during restart.

Literal Example
===============

This restarted-checkpoint block is illustrative. Its UUIDs and timestamps
are not real records.

.. code-block:: text

   ERF_PROVENANCE_BEGIN
   schema_version=1
   simulation_uuid=11111111-1111-4111-8111-111111111111
   execution_uuid=22222222-2222-4222-8222-222222222222
   parent_execution_uuid=11111111-1111-4111-8111-111111111111
   execution_lineage=11111111-1111-4111-8111-111111111111,22222222-2222-4222-8222-222222222222
   restart_generation=1
   lineage_complete=true
   lineage_status=complete
   source_checkpoint_uuid=aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa
   source_checkpoint_path=chk00010
   execution_start_utc=2026-07-11T19:00:00Z
   artifact_uuid=bbbbbbbb-bbbb-4bbb-8bbb-bbbbbbbbbbbb
   artifact_type=checkpoint
   artifact_step=20
   artifact_time_seconds=120
   artifact_created_utc=2026-07-11T19:05:00Z
   ERF_PROVENANCE_END

Legacy Checkpoints
==================

Older checkpoints have no provenance block. Missing ``job_info``, a missing
exact block, malformed fields, invalid UUIDs, failed lineage validation, an
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
``job_info`` file. ERF writes it to native checkpoints and primary native 2D
and 3D AMReX plotfiles. It does not write a sidecar.

ERF generates the execution UUID on the I/O rank and broadcasts the execution
record to all ranks. The I/O rank generates a separate artifact UUID when it
writes each ``job_info``. Artifact UUIDs are output metadata and are not
needed on non-I/O ranks. UUID generation does not use AMReX's scientific
random-number stream. UUIDs identify records but do not authenticate their
contents; they are not signatures or cryptographic evidence.

Scope and Limits
================

The checkpoint ``Header`` remains unchanged. NetCDF output, NetCDF column
files, 2D diagnostic JSON, subvolumes, staggered velocity plotfiles,
boundary-plane files, and particle-only files are outside this implementation.
This feature does not compute configuration hashes, digest input files or
executables, export W3C PROV or RO-Crate records, sign metadata, or reject a
restart based on provenance.
