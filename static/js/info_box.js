function raspCheckbox(target, textHtml, checked) {
    var checkbox = `
        <span id="radio-title" class="grey-text text-darken-1">RaSP Model<a
                    class="waves-effect waves-light btn modal-trigger grey lighten-2" href="#${target}"
                    title="Further information" style="float:right;">i</a></span>

        <div id="${target}" class="modal">
            <div class="modal-content">
                ${textHtml}
            </div>
            <div class="modal-footer">
                <a href="#!" class="modal-close waves-effect waves-green btn-flat">Close</a>
            </div>
        </div>

        <p>
            <label>
                <input type="checkbox" name="rasp-checkbox" ${checked} />
                <span>Calculate RaSP model (~10 min for each chain)</span>
            </label>
        </p>
    `;

    return checkbox;
}

function interfaceScoreCheckbox(target, textHtml, checked) {
    var checkbox = `
        <span id="radio-title" class="grey-text text-darken-1">Interface Score 
            <a class="waves-effect waves-light btn modal-trigger grey lighten-2" href="#${target}"
            title="Further information" style="float:right;">i</a>
        </span>

        <div id="${target}" class="modal">
            <div class="modal-content">
                ${textHtml}
            </div>
            <div class="modal-footer">
                <a href="#!" class="modal-close waves-effect waves-green btn-flat">Close</a>
            </div>
        </div>

        <p>
            <label>
                <input type="checkbox" name="ifscore_checkbox" ${checked} />
                <span>Calculate Interface score</span>
            </label>
        </p>
    `;

    return checkbox;
}

function minimizationCheckbox(target, htmlText) {
    var checkbox = `
        <span id="radio-title" class="grey-text text-darken-1">Minimization
            <a class="waves-effect waves-light btn modal-trigger grey lighten-2" href="#${target}"
                title="Further information" style="float:right;">i</a>
        </span>

        <div id="${target}" class="modal">
            <div class="modal-content">
                ${htmlText}
            </div>
            <div class="modal-footer">
                <a href="#!" class="modal-close waves-effect waves-green btn-flat">Close</a>
            </div>
        </div>

        <p>
            <label>
                <input type="radio" name="min-selector" value="long" checked required />
                <span>long (5-60 min, recommended)</span>
            </label>
        </p>
        <p>
            <label>
                <input type="radio" name="min-selector" value="short" />
                <span>quick (1-3 min)</span>
            </label>
        </p>
        <p>
            <label>
                <input type="radio" name="min-selector" value="none" />
                <span>none<span>
            </label>
        </p>
    `;

    return checkbox;
}