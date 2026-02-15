const CODON_TABLE = {
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
  TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G'
};

const analyzeBtn = document.getElementById('analyzeBtn');
const resultsSection = document.getElementById('results');
const errorEl = document.getElementById('error');

analyzeBtn.addEventListener('click', () => {
  const original = sanitizeSequence(document.getElementById('originalSequence').value);
  const mutated = sanitizeSequence(document.getElementById('mutatedSequence').value);

  if (!original || !mutated) {
    return showError('Please provide both original and mutated DNA sequences.');
  }
  if (!isValidDNA(original) || !isValidDNA(mutated)) {
    return showError('DNA can only contain A, T, C, and G characters.');
  }

  showError('');

  const alignment = alignSequences(original, mutated);
  const differences = summarizeDifferences(alignment);
  const mutationTypes = detectMutationTypes(differences);

  const originalProtein = translateDNA(original);
  const mutatedProtein = translateDNA(mutated);
  const proteinEffect = determineProteinEffect(originalProtein, mutatedProtein);

  const features = buildFeatures({ original, mutated, differences, proteinEffect, originalProtein, mutatedProtein });
  const prediction = predictImpact(features);

  renderResults({ differences, mutationTypes, originalProtein, mutatedProtein, proteinEffect, features, prediction });
});

function sanitizeSequence(sequence) {
  return sequence.toUpperCase().replace(/\s+/g, '');
}

function isValidDNA(sequence) {
  return /^[ATCG]+$/.test(sequence);
}

function showError(message) {
  errorEl.textContent = message;
}

function alignSequences(original, mutated) {
  const n = original.length;
  const m = mutated.length;
  const dp = Array.from({ length: n + 1 }, () => Array(m + 1).fill(0));
  const back = Array.from({ length: n + 1 }, () => Array(m + 1).fill(null));

  for (let i = 1; i <= n; i++) {
    dp[i][0] = i;
    back[i][0] = 'up';
  }
  for (let j = 1; j <= m; j++) {
    dp[0][j] = j;
    back[0][j] = 'left';
  }

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const matchOrSubCost = dp[i - 1][j - 1] + (original[i - 1] === mutated[j - 1] ? 0 : 1);
      const deleteCost = dp[i - 1][j] + 1;
      const insertCost = dp[i][j - 1] + 1;

      const best = Math.min(matchOrSubCost, deleteCost, insertCost);
      dp[i][j] = best;

      if (best === matchOrSubCost) back[i][j] = 'diag';
      else if (best === deleteCost) back[i][j] = 'up';
      else back[i][j] = 'left';
    }
  }

  const operations = [];
  let i = n;
  let j = m;
  while (i > 0 || j > 0) {
    const move = back[i][j];
    if (move === 'diag') {
      operations.push({
        type: original[i - 1] === mutated[j - 1] ? 'match' : 'substitution',
        from: original[i - 1],
        to: mutated[j - 1],
        originalIndex: i,
        mutatedIndex: j
      });
      i--;
      j--;
    } else if (move === 'up') {
      operations.push({
        type: 'deletion',
        from: original[i - 1],
        to: '-',
        originalIndex: i,
        mutatedIndex: j
      });
      i--;
    } else {
      operations.push({
        type: 'insertion',
        from: '-',
        to: mutated[j - 1],
        originalIndex: i,
        mutatedIndex: j
      });
      j--;
    }
  }

  operations.reverse();
  return operations;
}

function summarizeDifferences(operations) {
  return operations.filter((op) => op.type !== 'match');
}

function detectMutationTypes(differences) {
  const found = new Set(differences.map((item) => item.type));
  if (found.size === 0) return ['No mutation detected'];
  return Array.from(found).map((type) => type[0].toUpperCase() + type.slice(1));
}

function translateDNA(dna) {
  let protein = '';
  for (let i = 0; i + 2 < dna.length; i += 3) {
    const codon = dna.slice(i, i + 3);
    protein += CODON_TABLE[codon] || 'X';
  }
  return protein;
}

function determineProteinEffect(originalProtein, mutatedProtein) {
  if (originalProtein === mutatedProtein) return 'Silent mutation';

  const stopOriginal = originalProtein.indexOf('*');
  const stopMutated = mutatedProtein.indexOf('*');

  const hasPrematureStop =
    stopMutated !== -1 && (stopOriginal === -1 || stopMutated < stopOriginal);

  if (hasPrematureStop) return 'Nonsense mutation';

  return 'Missense mutation';
}

function buildFeatures({ original, mutated, differences, proteinEffect, originalProtein, mutatedProtein }) {
  const substitutions = differences.filter((d) => d.type === 'substitution').length;
  const insertions = differences.filter((d) => d.type === 'insertion').length;
  const deletions = differences.filter((d) => d.type === 'deletion').length;
  const indels = insertions + deletions;
  const frameshift = (mutated.length - original.length) % 3 !== 0 ? 1 : 0;

  const comparableLength = Math.max(1, Math.min(originalProtein.length, mutatedProtein.length));
  let proteinChanges = 0;
  for (let i = 0; i < comparableLength; i++) {
    if (originalProtein[i] !== mutatedProtein[i]) proteinChanges++;
  }
  const proteinChangeRatio = proteinChanges / comparableLength;

  return {
    substitutions,
    insertions,
    deletions,
    indels,
    frameshift,
    silent: proteinEffect === 'Silent mutation' ? 1 : 0,
    missense: proteinEffect === 'Missense mutation' ? 1 : 0,
    nonsense: proteinEffect === 'Nonsense mutation' ? 1 : 0,
    proteinChangeRatio
  };
}

function predictImpact(features) {
  const linearScore =
    -1.4 +
    0.35 * features.substitutions +
    0.5 * features.indels +
    1.2 * features.frameshift +
    1.4 * features.nonsense +
    0.5 * features.missense -
    0.9 * features.silent +
    1.6 * features.proteinChangeRatio;

  const probability = 1 / (1 + Math.exp(-linearScore));
  const label = probability >= 0.55 ? 'Potentially Harmful' : 'Likely Neutral';

  return { probability, label };
}

function renderResults({ differences, mutationTypes, originalProtein, mutatedProtein, proteinEffect, features, prediction }) {
  resultsSection.hidden = false;

  document.getElementById('differenceCount').textContent = String(differences.length);
  document.getElementById('mutationTypes').textContent = mutationTypes.join(', ');
  document.getElementById('originalProtein').textContent = originalProtein || '(no complete codon translated)';
  document.getElementById('mutatedProtein').textContent = mutatedProtein || '(no complete codon translated)';
  document.getElementById('proteinEffect').textContent = proteinEffect;

  const differenceList = document.getElementById('differenceList');
  differenceList.innerHTML = '';
  if (differences.length === 0) {
    differenceList.innerHTML = '<li>No difference found.</li>';
  } else {
    for (const diff of differences) {
      const item = document.createElement('li');
      item.textContent = `${diff.type.toUpperCase()} at Original #${diff.originalIndex}, Mutated #${diff.mutatedIndex}: ${diff.from} → ${diff.to}`;
      differenceList.appendChild(item);
    }
  }

  document.getElementById('predictionLabel').textContent = prediction.label;
  document.getElementById('predictionScore').textContent = `${(prediction.probability * 100).toFixed(1)}%`;

  const featureList = document.getElementById('featureList');
  featureList.innerHTML = '';
  for (const [key, value] of Object.entries(features)) {
    const li = document.createElement('li');
    li.textContent = `${key}: ${typeof value === 'number' ? value.toFixed(3).replace(/\.000$/, '') : value}`;
    featureList.appendChild(li);
  }
}
