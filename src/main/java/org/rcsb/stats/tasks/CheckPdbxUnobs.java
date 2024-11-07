package org.rcsb.stats.tasks;

import org.rcsb.cif.CifIO;
import org.rcsb.cif.ParsingException;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.*;
import org.rcsb.stats.Constants;
import org.rcsb.stats.Helpers;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Obtain the number of non-hydrogen atoms as described by the `atom_site` content of an mmCIF file.
 */
public class CheckPdbxUnobs {
    static final Logger logger = LoggerFactory.getLogger(CheckPdbxUnobs.class);
    final AtomicInteger counter = new AtomicInteger();

    public static void main(String[] args) throws IOException {
        CheckPdbxUnobs me = new CheckPdbxUnobs();
        me.countHeavyAtoms(me.fetchStructureData("5MG3"));
        new CheckPdbxUnobs().computeStats();
    }

    void computeStats() throws IOException {
        // obtain set of all known identifiers (1ABC, 1ABD, ...)
        Set<String> identifiers = Helpers.getAllIdentifiers(Constants.ResultsContentType.EXPERIMENTAL);

        // traverse all identifiers
        long heavyAtomCount = identifiers.parallelStream().peek(this::logProgress)
                // load the {@link CifFile} for each identifier
                .map(this::fetchStructureData)
                // process each structure: count the number of non-hydrogen atoms
                .mapToLong(this::countHeavyAtoms)
                // aggregate as sum
                .sum();

        logger.info("There are {} heavy (non-hydrogen) atoms in {} PDB structures", Helpers.formatNumber(heavyAtomCount), Helpers.formatNumber(counter.get()));
    }

    /**
     * Obtain structure for a single identifier from models.rcsb.org.
     * @param identifier datum to load
     * @return structure data in the form of a {@link CifFile}
     */
    CifFile fetchStructureData(String identifier) {
        try {
            URL url = new URL(String.format(Constants.BCIF_SOURCE, identifier));
            // other CifIO methods allow reading from Paths or byte streams -- methods for writing can be found there too
            return CifIO.readFromURL(url);
        } catch (IOException e) {
            logger.warn("Failed to pull structure data for {}", identifier);
            throw new UncheckedIOException(e);
        } catch (ParsingException e) {
            logger.warn("Failed to parse structure data for {}", identifier);
            throw e;
        }
    }

    /**
     * Process a CIF file (obtained from @link{CifIO#readFromUrl()}).
     * @param cifFile source data
     * @return the count of non-hydrogen atoms
     */
    long countHeavyAtoms(CifFile cifFile) {
        Map<String, Integer> entityIdToLength = new HashMap<>();
        Map<String, String> asymIdsToEntId  = new HashMap<>();
        MmCifBlock block = cifFile
                // optional: apply mmCIF schema to get schema definitions and types
                .as(StandardSchemata.MMCIF)
                // CIF files may have multiple blocks of data, the PDB archive only makes use of the 1st
                .getFirstBlock();

        String entryId = block.getStruct().getEntryId().values().findFirst().orElse("");
        boolean isDiffraction = block.getExptl().getMethod().values().anyMatch(m -> m.equals("X-RAY DIFFRACTION") || m.equals("NEUTRON DIFFRACTION"));

        PdbxNmrEnsemble ensemble = block.getPdbxNmrEnsemble();
        int numModels = ensemble.getConformersSubmittedTotalNumber().values().findAny().orElse(1);

        EntityPolySeq entityPoly = block.getEntityPolySeq();
        for (int rowIndex = 0; rowIndex < entityPoly.getRowCount(); rowIndex++) {
            String entId = entityPoly.getEntityId().get(rowIndex);
            // this essentially achieves getting the last value of an entity. A bit heavy to put for every row but an easier implementation ;)
            // the alternative of counting all rows doesn't work because of microhetereogeneity (e.g. 3ZE8)
            entityIdToLength.put(entId, entityPoly.getNum().get(rowIndex));
        }

        AtomSite atomSite = block.getAtomSite();
        Map<String, TreeSet<Integer>> asymIdToResNums = new HashMap<>();
        double totalResOcc = 0.0;
        int prevSeqId = -1;
        String prevAsymId = null;
        double minOcc = -1.0; // this equates to ignoring occupancy for non-diffraction methods
        if (isDiffraction) minOcc = 0.00001;
        for (int rowIndex = 0; rowIndex < atomSite.getRowCount(); rowIndex++) {
            String entId = atomSite.getLabelEntityId().get(rowIndex);
            String asymId = atomSite.getLabelAsymId().get(rowIndex);
            asymIdsToEntId.putIfAbsent(asymId, entId);
            if (!entityIdToLength.containsKey(entId)) continue; // not a polymer
            // init to empty sets: if it happens that the loop through atom site doesn't find anything for an entity, then this will reflect it (e.g. 5MG3)
            asymIdToResNums.putIfAbsent(asymId, new TreeSet<>());
            int seqId = atomSite.getLabelSeqId().get(rowIndex);
            double occ = atomSite.getOccupancy().get(rowIndex);

            if (prevSeqId>0 && (seqId!=prevSeqId || !asymId.equals(prevAsymId)) && totalResOcc > minOcc) {
                asymIdToResNums.get(prevAsymId).add(prevSeqId);
            }

            if (prevSeqId>0 && seqId!=prevSeqId) {
                totalResOcc = 0.0;
            }
            totalResOcc += occ;

            prevSeqId = seqId;
            prevAsymId = asymId;
        }
        if (totalResOcc > minOcc) {
            asymIdToResNums.get(prevAsymId).add(prevSeqId);
        }

        Map<String, List<String>> entIdToAsymIds = new HashMap<>();
        asymIdsToEntId.forEach((key, value) -> entIdToAsymIds.computeIfAbsent(value, k -> new ArrayList<>()).add(key));

        Map<String, Integer> asymIdToUnmodeledCount = new HashMap<>();
        for (String entId : entityIdToLength.keySet()) {
            for (String asymId : entIdToAsymIds.get(entId)) {
                for (int i = 1; i <= entityIdToLength.get(entId); i++) {
                    if (!asymIdToResNums.get(asymId).contains(i)) {
                        asymIdToUnmodeledCount.merge(asymId, 1, Integer::sum);
                    }
                }
            }
        }

        PdbxUnobsOrZeroOccResidues pdbxUnobsRes = block.getPdbxUnobsOrZeroOccResidues();
        Map<String, Integer> asymIdToPdbUnobsCount = new HashMap<>();
        for (int rowIndex = 0; rowIndex < pdbxUnobsRes.getRowCount(); rowIndex++) {
            String asymId = pdbxUnobsRes.getLabelAsymId().get(rowIndex);
            asymIdToPdbUnobsCount.merge(asymId, 1, Integer::sum);
        }

        for (String entId : entityIdToLength.keySet()) {
            for (String asymId : entIdToAsymIds.get(entId)) {
                int unmodeledCount = Optional.ofNullable(asymIdToUnmodeledCount.get(asymId)).orElse(0) * numModels;
                int unobsCount = Optional.ofNullable(asymIdToPdbUnobsCount.get(asymId)).orElse(0);
                if (unmodeledCount != unobsCount) {
                    String multiModMsg = "";
                    if (numModels>1) multiModMsg = "Note entry has " + numModels + " models. ";
                    if (!isDiffraction) multiModMsg += "Not a diffraction entry.";
                    logger.info("Different count for entry {}, asym {}, unmodeled {}, unobs {}. {}", entryId, asymId, unmodeledCount, unobsCount, multiModMsg);
                }
            }
        }

        return 0;
    }

    final Set<String> HYDROGEN_ATOMS = Set.of("H", "D", "T");

    /**
     * Filters for non-hydrogen atoms based on their `atom_site.type_symbol`.
     * @param typeSymbol element of this atom
     * @return false if this is hydrogen
     */
    boolean isHeavyAtom(String typeSymbol) {
        return !HYDROGEN_ATOMS.contains(typeSymbol);
    }

    /**
     * Indicate that 10,000 entries have been processed.
     */
    void logProgress(Object ignored) {
        if (counter.incrementAndGet() % 10000 == 0) {
            logger.info("Processed {} entries", Helpers.formatNumber(counter.get()));
        }
    }
}
