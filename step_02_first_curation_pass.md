# Step 2: First Curation Pass

Now that the operational definitions are fixed, the next task is not sequence
collection yet.

The next task is to produce the **first 10 to 20 core rows** in the curated
taxon table.

## Immediate workflow

1. Pick one class first.
2. Screen a small number of trusted sources.
3. Enter every candidate into the screening log.
4. Transfer only supported rows into the core taxon table.
5. Mark each row as `yes`, `review`, or `no`.

## Recommended order

Start with `thermophiles` first.

Why this is a good first pass:

- OGT values are often reported explicitly
- the threshold is simple
- `hyperthermophile` is a clean subclass
- the literature language is usually more standardized

After that, do the same for `halophiles`.

## Small, realistic target

For the first pass, aim for:

- 10 thermophile rows
- 10 halophile rows

That is enough to stress-test:

- your schema
- your evidence vocabulary
- your inclusion rules
- your notes style

## Decision rule reminder

- `yes`: threshold satisfied and evidence is `experimental`,
  `curated_database`, or `review`
- `review`: only `genus_inference`, conflicting values, ambiguous identity, or
  values very close to the threshold
- `no`: threshold not satisfied, contradiction exists, or only
  `biome_inference` is available

## What not to do yet

Do not do these yet:

- collect raw protein FASTA files
- merge AFESM candidates into the core set
- add more extremophile classes
- worry about final model splits

## Output of this step

At the end of this step, you should have:

- a filled screening log
- a first curated taxon table
- a better sense of which sources are most efficient

Only then should you move to sequence retrieval.
