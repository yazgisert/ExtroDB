# DEPRECATED — this file has been renamed and fixed.
# Use prepare_combined_fasta.py instead.
#
# Changes from this version:
#   - Renamed (removed step1_ prefix)
#   - Added docstring explaining Split 1 collapsing (happens upstream in
#     build_searchable_extremophile_index.py, not here)
#   - Fixed halophile subclass matching in derive_split2_label:
#       "slight halophile"    -> "halophile (slight)"
#       "moderate halophile"  -> "halophile (moderate)"
#       "extreme halophile"   -> "halophile (extreme)"
#   - Added runtime assertion that every selected organism has a valid Split 1 label
#   - derive_split2_label now raises ValueError instead of silently returning "unknown"
raise SystemExit("Deprecated. Run prepare_combined_fasta.py instead.")
