async function fetchClustalWFile(clustalWFilePath) {
    return fetch(clustalWFilePath)
        .then(response => response.text())
        .then(data => {
            return data;
        })
        .catch(error => {
            console.error('Error fetching ClustalW file:', error);
            throw error;
        });
}

function extractChainIdsClustal(clustalWData) {
    const chainIdRegex = /^(?!CLUSTAL)(\S+)/gm;
    const uniqueChainIds = new Set();
    let match;

    while ((match = chainIdRegex.exec(clustalWData)) !== null) {
        uniqueChainIds.add(match[1]);
    }
    const chainIds = Array.from(uniqueChainIds);

    return chainIds;
}