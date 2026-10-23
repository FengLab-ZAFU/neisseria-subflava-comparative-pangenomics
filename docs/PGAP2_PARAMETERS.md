# PGAP2 parameters

The main pangenome used PGAP2 v2.4.1 with `dup_id=0.99` and `orth_id=0.98`. Hierarchical clustering used a coverage threshold of 0.98 for both sequences. Subsequent protein-alignment filters used `AS=0.60`, `AL=0.60`, `LD=0.60` and `para_id=0.70`, with an E-value threshold of 1e-5 and up to 2,000 targets per query.

Clustering coverage and post-alignment coverage are separate parameters. AS and AL use alignment length divided by the shorter and longer sequence lengths, respectively. LD is the shorter/longer sequence-length ratio. In the installed PGAP2 v2.4.1 source, clustering coverage was fixed at 0.98 and was not exposed as a command-line option.

The values were checked against the main run's saved parameter record and installed source. The original `partition.log` and installed source are not bundled. The sensitivity run changed `orth_id` to 0.95; its full saved parameter record was not checked again during release preparation.
