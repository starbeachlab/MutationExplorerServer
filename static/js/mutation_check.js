function getChainRanges(rawRanges) {
    var chains = rawRanges.replace(/\s/g, "").replace(/[,]$/, "").split(",");

    var ranges = [];
    for (c of chains) {
        r = c.split(/[:-]/);

        // if it is a negative range, the string will split into 4
        if (r.length == 4 && r[1] == '') {
            r.splice(1, 1);
            r[1] = `-${r[1]}`;
        }
        ranges.push(r);
    }

    return ranges
}

function mutationSyntaxCorrect(inputId) {
    muts = $(`#${inputId}`).val().replace(/ /g, "");
    if (!muts) {
        return true;
    }
    return /^([A-Z]:-?[0-9]+[ACDEFGHIKLMNPQRSTVWY])(,([A-Z]:-?[0-9]+[ACDEFGHIKLMNPQRSTVWY]))*$/.test(muts)
}

function checkMutation(mutation, ranges) {
    var chain = mutation[0];
    var position = parseInt(mutation.match(/-?\d+(\.\d+)?/));
    for (r of ranges) {
        if (r[0] == chain) {
            var start = parseInt(r[1]);
            var end = parseInt(r[2]);
            if (start <= position && position <= end) {
                return true;
            }
        }
    }
    return false;
}

function mutationChainCheck(inputId, rawRanges) {
    var mutations = $(`#${inputId}`).val().replace(/\s/g, "");
    if (mutations.length === 0) {
        return true;
    }

    mutations = mutations.split(",");
    ranges = getChainRanges(rawRanges)

    for (m of mutations) {
        if (!checkMutation(m, ranges)) {
            return false
        }
    }

    return true;
};

function mutationCheck(inputField, errorField, rawRanges, e) {
    if (!mutationSyntaxCorrect(inputField)) {
        e.preventDefault();
        console.log("incorrect mutation syntax");
        $(`#${inputField}`).addClass("invalid");
        $(`#${errorField}`).html("Invalid syntax. If your syntax is correct, please use upper case letters. See info button for more information");
    } else if (!mutationChainCheck(inputField, rawRanges)) {
        e.preventDefault();
        console.log("incorrect mutation position (chain or residue)");
        $(`#${inputField}`).addClass("invalid");
        $(`#${errorField}`).html("Invalid chain or residue position. See info button for more information");
    } 
}