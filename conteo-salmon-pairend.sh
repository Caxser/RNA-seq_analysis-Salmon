#!/bin/bash

# Ruta al índice de Salmon
INDICE="/mnt/d/Proyecto-maestria/indice"

# Directorio con los archivos fastq
FASTQ_DIR="/mnt/d/Proyecto-maestria/fastq_clean"

# Directorio de salida
SALIDA_DIR="/mnt/d/Proyecto-maestria/conteos-salmon"

# Número de hilos
HILOS=10

# Crear carpeta de salida si no existe
mkdir -p "$SALIDA_DIR"

# Recorrer todos los archivos *_1.fastq.gz
for R1 in "$FASTQ_DIR"/*/*_R1.clean.fastq.gz; do
    # Extraer nombre base (por ejemplo, SRR12038075)
    BASE=$(basename "$R1" _R1.fastq.gz)
    BASE2=${BASE%%_*}

    DIREC=$(dirname "$R1")   
    R2="${DIREC}/${BASE2}_R2.clean.fastq.gz"
    echo "Queda asi:  $R2"
    
    # Validar que ambos archivos existen
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "Procesando muestra: $BASE"
        
        salmon quant -i "$INDICE" \
            -l A \
            -1 "$R1" \
            -2 "$R2" \
            -p "$HILOS" \
            --validateMappings \
            --seqBias --gcBias --posBias \
            -o "${SALIDA_DIR}/${BASE2}"
    else
        echo "⚠️  Archivos faltantes para la muestra: $BASE"
    fi
done
