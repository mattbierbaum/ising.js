// all of the global variables for dynamics
var gpx_black = null;
var gpx_white = null;
var gpx_size = 0;
var canvasN = 512;
var gbuffer;
var gbufferdata;

var gboard = null;
var gN = 256;
var gT = 2.26918531421;
var gfield = 0;

var ge_avg, ge_var, gm_avg, gm_var;
var gtable_de;
var gtable_doflip;
var gtable_flipprob;
var wolfp = 1-Math.exp(-2./gT);

var gt = 0;
var times = [];
var gtimeseries_energy = [];
var gtimeseries_mag = [];
var gtimeseries_eavg = [];
var gtimeseries_mavg = [];
var genergy = 0;
var gmag = 0;
var frame = 0;
var sweeps = 0;
var lasttime = 0;

// display variables
var c, c2;
var ctx;
var ctxgraph;
var empty;
var frameskip = 1;
var onefill = 0;
var dodraw = true;
var gh = 150;
var gw = 370;

function rgb(r,g,b) {
    return 'rgb('+r+','+g+','+b+')';
}
function log10(val) {
    return Math.log(val) / Math.LN10;
}

function toFixed(value, precision, negspace) {
    negspace = typeof negspace !== 'undefined' ? negspace : '';
    var precision = precision || 0;
    var sneg = (value < 0) ? "-" : negspace;
    var neg = value < 0;
    var power = Math.pow(10, precision);
    var value = Math.round(value * power);
    var integral = String(Math.abs((neg ? Math.ceil : Math.floor)(value/power)));
    var fraction = String((neg ? -value : value) % power);
    var padding = new Array(Math.max(precision - fraction.length, 0) + 1).join('0');
    return sneg + (precision ? integral + '.' +  padding + fraction : integral);
}

function init_board(N, board){
    gt = 0;
    gboard = [];
    gN = N;

    if (board !== null){
        for (var i=0; i<gN*gN; i++)
            gboard[i] = board[i]; 
    } else {
        for (var i=0; i<gN*gN; i++)
            gboard[i] = 2*Math.floor(Math.random()*2) - 1;
    }

    gpx_size = canvasN/gN;
    display_board(gN, gboard);
    draw_all();

    init_measurements();
}

function put_pixel(x, y, size, color){
    var xoff = x*size;
    var yoff = y*size;
    for (var i=0; i<size; i++){
        for (var j=0; j<size; j++){
            var ind = ((yoff+j)*gN*size + xoff+ i)*4;
            var c = (color+1)/2 * 255;
            gbufferdata[ind+0] = c;
            gbufferdata[ind+1] = c;
            gbufferdata[ind+2] = c;
            gbufferdata[ind+3] = 255;
        }
    }
}

function display_board(N, board){
    for (var i=0; i<N; i++){
        for (var j=0; j<N; j++){
            put_pixel(i, j, gpx_size, board[i+j*N]);             
        }
    }
}

function energy(x, y, N, b){
    return -b[x+y*N]*(b[x + ((y+1).mod(N))*N] + 
        b[x + ((y-1).mod(N))*N] + 
        b[(x+1).mod(N) + y*N] + 
        b[(x-1).mod(N) + y*N] + gfield);
}


function energy_difference(x, y, N, b){
    return -2*energy(x,y,N,b);
}

function neighborCount(x, y, N, b){
    var i = x+y*N;
    return 1*(b[x + ((y+1).mod(N))*N] > 0) + 
        1*(b[x + ((y-1).mod(N))*N] > 0) + 
        1*(b[(x+1).mod(N) + y*N] > 0) +
        1*(b[(x-1).mod(N) + y*N] > 0);
}

function update_metropolis(){
    var x = Math.floor(Math.random()*gN);
    var y = Math.floor(Math.random()*gN);
    var ind = x + y*gN;
    var de = energy_difference(x, y, gN, gboard);
    if (de <= 0 || Math.random() < Math.exp(-de / gT)){
        gboard[ind] = -gboard[ind];

        if (!onefill)
            put_pixel(x, y, gpx_size, gboard[x+y*gN]);

        genergy += 1.0*de/(gN*gN);
        gmag += 2.0*gboard[ind]/(gN*gN);
    }
    gt += 1.0/(gN*gN);
}

function update_wolff() {
    // pick random site to start
    var x = Math.floor(Math.random()*gN);
    var y = Math.floor(Math.random()*gN);
    var ind = x + y*gN;

    // get the initial state and seed
    // the cluster
    var sites = [ind];
    var state = gboard[ind];

    var is_good = function (next_ind) {
        next_ind = Number(next_ind);
        return (
                (!(next_ind in cluster)) &&
                (Math.random() < wolfp) &&
                (gboard[next_ind]==state)
               )
    }

    var cluster = {};
    var frontier = {};
    cluster[ind] = 1;
    frontier[ind] = 1;
    var newfrontier = {};
    var next_ind = 0;

    while (Object.keys(frontier).length > 0) {
        newfrontier = {};

        for (var current_ind in frontier) {
            current_ind = Number(current_ind);
            x = current_ind.mod(gN);
            y = Math.floor( current_ind / gN );

            // do each neighbor
            next_ind = x + ((y+1).mod(gN))*gN;
            if (is_good(next_ind)) {
                newfrontier[next_ind] = 1;
                cluster[next_ind] = 1;
            }
            next_ind = x + ((y-1).mod(gN))*gN;
            if (is_good(next_ind)) {
                newfrontier[next_ind] = 1;
                cluster[next_ind] = 1;
            }
            next_ind = (x+1).mod(gN) + y*gN;
            if (is_good(next_ind)) {
                newfrontier[next_ind] = 1;
                cluster[next_ind] = 1;
            }
            next_ind = (x-1).mod(gN) + y*gN;
            if (is_good(next_ind)) {
                newfrontier[next_ind] = 1;
                cluster[next_ind] = 1;
            }
        }
        frontier = newfrontier;
    }
    // having built the cluster, determine the probability of flipping
    var clussize = Object.keys(cluster).length;
    var ds = 2 * state * clussize * gfield;

    if ( (ds <= 0) || (Math.random() < Math.exp(-ds)) ) {
        // flip the cluster
        for (var ind in cluster) {
            ind = Number(ind);
            x = ind % gN;
            y = Math.floor( ind / gN );
            gboard[ind] = -state;
            put_pixel(x,y, gpx_size, -state);

            de = energy_difference(x, y, gN, gboard);
            genergy += 1.0*de/(gN*gN);
            gmag += 2.0*gboard[ind]/(gN*gN);
        }
    }
    gt += clussize / (gN*gN);
}


var update_func = "metropolis";
function update() {
    if (update_func=="metropolis") {
        update_metropolis();
    }
    else if (update_func=="wolff") {
        update_wolff();
    }
}

function push_measurement(t, e, m){
    times.push(t);
    gtimeseries_energy.push(e);
    gtimeseries_mag.push(m);

    n = times.length;
    ge0 = ge_avg;
    gm0 = gm_avg;

    // welford's algorithm
    ge_avg = ge_avg + (e - ge_avg)/n;
    ge_var = ((n-1)*ge_var + (e - ge_avg)*(e - ge0)) / n;
    gm_avg = gm_avg + (m - gm_avg)/n;
    gm_var = ((n-1)*gm_var + (m - gm_avg)*(m - gm0)) / n;
    gtimeseries_eavg.push(ge_avg);
    gtimeseries_mavg.push(gm_avg);

    sps = 1000.0*sweeps/(Date.now() - lasttime);
}

function init_measurements(){
    frame = 0;
    sweeps = 0;
    gt = 0;
    ge_avg = ge_var = gm_avg = gm_var = 0;
    times = [];
    gtimeseries_energy = [];
    gtimeseries_mag = [];
    gtimeseries_eavg = [];
    gtimeseries_mavg = [];
    lasttime = Date.now();
    reset_measurements();
    push_measurement(gt, genergy, gmag);
}

function reset_measurements(){
    genergy = 0;
    gmag = 0;
    ge_avg = ge_var = gm_avg = gm_var = 0;

    for (var i=0; i<gN; i++){
    for (var j=0; j<gN; j++){
        genergy += energy(i,j,gN, gboard);
        gmag += gboard[i+gN*j];
    }}
    genergy /= gN*gN*2;
    gmag /= gN*gN;

    ge_avg = genergy;
    gm_avg = gmag;
    update_measurements_labels();
}

function update_measurements_labels(){
    lblt = document.getElementById('label_time');
    lble = document.getElementById('label_energy');
    lblm = document.getElementById('label_mag');

    lblt.innerHTML = "time = "+toFixed(gt, 4, ' ')+"   sweeps/sec = "+toFixed(sps, 3, ' ');
    lble.innerHTML = "e = "+toFixed(genergy, 5, ' ');
    lblm.innerHTML = "m = "+toFixed(gmag, 5, ' ');

    lble.innerHTML += "   &lt;e&gt; = "+toFixed(ge_avg, 5, ' ');
    lblm.innerHTML += "   &lt;m&gt; = "+toFixed(gm_avg, 5, ' ');

    lble.innerHTML += "   Var(e) = "+toFixed(ge_var, 7, ' ');
    lblm.innerHTML += "   Var(m) = "+toFixed(gm_var, 7, ' ');
}

function hidden_link_download(uri, filename){
    var link = document.createElement('a');
    link.href = uri;
    link.style.display = 'none';
    link.download = filename
    link.id = 'templink';
    document.body.appendChild(link);
    document.getElementById('templink').click();
    document.body.removeChild(document.getElementById('templink'));
}

function download_measurements(){
    var csv = "data:text/csv;charset=utf-8,";
    csv += "# time, energy per spin, magnetization per spin\n";
    for (var i=0; i<times.length; i++){
        csv += times[i]+", ";
        csv += gtimeseries_energy[i]+", ";
        csv += gtimeseries_mag[i]+"\n";
    }
    var encoded = encodeURI(csv);
    hidden_link_download(encoded, 'ising-data.txt');
}

function download_field(){
    uri = c.toDataURL("image/png");
    hidden_link_download(uri, 'ising-field.png');
}

function download_graph(){
    uri = c2.toDataURL("image/png");
    hidden_link_download(uri, 'ising-graph.png');
}

function draw_all(){
    if (onefill)
        display_board(gN, gboard);

    gbuffer.data = gbufferdata;
    ctx.putImageData(gbuffer, 0, 0);
    push_measurement(gt, genergy, gmag);
    update_measurements_labels();
    draw_graph();
}

function draw_graph(){
    cleargraph();
    if (graph_type != "none"){
        x = times;
        if (graph_type == "energy")
            y = gtimeseries_energy;
        if (graph_type == "mag")
            y = gtimeseries_mag;
        if (graph_type == "eavg")
            y = gtimeseries_eavg;
        if (graph_type == "mavg")
            y = gtimeseries_mavg;

        draw_series_graph(x, y);
    }
}


/*======================================================================
  the javascript interface stuff
=========================================================================*/
function dotextbox(id){
    idt = id+"_input";
    document.getElementById(id).style.display = 'none';
    document.getElementById(idt).style.display = 'inline';
    document.getElementById(idt).value = document.getElementById(id).innerHTML;
    document.getElementById(idt).focus();
}

function undotextbox(id){
    idt = id.replace("_input", "");
    document.getElementById(idt).style.display = 'inline';
    document.getElementById(id).style.display = 'none';
}

function update_temp(){
    min = document.getElementById('temp').min;
    gTval = parseFloat(document.getElementById('temp').value);
    if (gTval <= min)
        gT = 0;
    else
        gT = Math.pow(10, gTval);
    document.getElementById('label_temp').innerHTML = toFixed(gT,6);
    calculateFlipTable(gT);
}
function update_field(){
    gfield = parseFloat(document.getElementById('field').value);
    document.getElementById('label_field').innerHTML = toFixed(gfield,6);
    calculateFlipTable(gT);
    reset_measurements();
}
function update_frames(){
    frameval = parseFloat(document.getElementById('frames').value);
    if (update_func=='metropolis') {
        frameskip = Math.pow(10, frameval);
        onefill = frameskip > 2*gN*gN ? 1 : 0;
    } else {
      frameskip = frameval;
        onefill = frameskip > 3 ? 1 : 0;
    }
    document.getElementById('label_frames').innerHTML = toFixed(frameskip,6);
}

function update_display(){
    document.getElementById('label_temp').innerHTML = toFixed(gT,6);
    document.getElementById('label_field').innerHTML = toFixed(gfield,6);
    document.getElementById('label_frames').innerHTML = toFixed(frameskip,6);
}

function update_method() {
    var frame_slider = document.getElementById('frames');
    var frame_label = document.getElementById('label_frames');
    if (document.getElementById('method_wolff').checked) {
        update_func = 'wolff';
        frameskip = 2;
        frame_label.innerHTML = toFixed(frameskip,0);
        frame_slider.step = 1;
        frame_slider.max=20;
        frame_slider.min=1;
        frame_slider.value = frameskip;
    } else  {
        update_func = 'metropolis';
        frameskip = Math.pow(10.,0);
        frame_label.innerHTML = toFixed(1,6);
        frame_slider.step = 0.01;
        frame_slider.max=2;
        frame_slider.min=-2;
        frame_slider.value = 0;
    }
}

function update_pause(){
    if (dodraw == true){
        document.getElementById('pause').value = 'Start';
        dodraw = false;
    } else {
        document.getElementById('pause').value = 'Pause';
        requestAnimationFrame(tick, c);
        dodraw = true;
    }
}

function update_restart(){
    init_board(gN, null);
}

function update_step(){
    if (dodraw)
        update_pause();

    if (update_func == 'metropolis'){
        for (var i=0; i<gN*gN; i++)
            update();
    } else {
        for (var i=0; i<Math.sqrt(gN); i++)
            update();
    }
    draw_all();
}

/*===============================================================================
 * graphing
 *=============================================================================*/
var xaxis = 40;
var yaxis = 5;
function x2px(x, xmin, dx) {return ((x - xmin) / dx) * (gw - xaxis) + xaxis; }
function y2px(y, ymin, dy) {return gh - ((y - ymin) / dy * (gh - 2*yaxis) + yaxis); }
var graph_type = "energy";

function draw_series_graph(xl, yl){
    var xllength = xl.length;
    var xmax, xmin, ymax, ymin;
    xmax = ymax = -1e10; xmin = ymin = 1e10;
    var skip = 1;//Math.floor(1+(xllength/gw));

    for (var i=0; i<xllength; i+=skip){
        if (xl[i] < xmin) xmin = xl[i];
        if (xl[i] > xmax) xmax = xl[i];
        if (yl[i] < ymin) ymin = yl[i];
        if (yl[i] > ymax) ymax = yl[i];
    }

    var dx = xmax - xmin;
    var dy = ymax - ymin;

    var oom_x = Math.abs(dx)<1?Math.round(log10(dx)):Math.floor(log10(dx));
    var oom_y = Math.abs(dy)<1?Math.round(log10(dy)):Math.floor(log10(dy));
    var pow10_x = Math.pow(10, oom_x);
    var pow10_y = Math.pow(10, oom_y);

    var idx = Math.floor(dx / pow10_x);
    var idy = Math.floor(dy / pow10_y);

    if (idx < 1) idx = 10; 
    if (idy < 1) idy = 10; 

    if (idx == 1 || idx == 2)
        idx *= 5;
    if (idy == 1 || idy == 2)
        idy *= 5;

    xtic_major = dx/idx;
    ytic_major = dy/idy;
    xtic_minor = xtic_major/5;
    ytic_minor = ytic_major/5;

    /*ymin = Math.floor(ymin/ytic_major)*ytic_major;
    ymax = Math.ceil(ymax/ytic_major)*ytic_major;*/

    ctxgraph.font='12px sans-serif';
    ctxgraph.fillStyle='rgba(0,0,0,1)';

    ctxgraph.beginPath();
    ctxgraph.moveTo(xaxis, 0);
    ctxgraph.lineTo(xaxis, gh);
    ctxgraph.stroke();

    ctxgraph.beginPath();
    ctxgraph.moveTo(xaxis, y2px(0, ymin, dy));
    ctxgraph.lineTo(gw, y2px(0, ymin, dy));
    ctxgraph.stroke();

    for (var i=-idy; i<=idy; i++){
        y = y2px(i*ytic_major+(ymin+ymax)/2, ymin, dy);
        ctxgraph.beginPath();
        ctxgraph.moveTo(xaxis-5, y);
        ctxgraph.lineTo(xaxis, y);
        ctxgraph.stroke();
        ctxgraph.fillText(toFixed(i*ytic_major+(ymin+ymax)/2, 3), 0, y+4);
    }

    for (var i=-idy*5; i<=idy*5; i++){
        y = y2px(i*ytic_minor+(ymin+ymax)/2, ymin, dy);
        ctxgraph.beginPath();
        ctxgraph.moveTo(xaxis-2, y);
        ctxgraph.lineTo(xaxis, y);
        ctxgraph.stroke();
    }

    for (var i=0; i<xllength-skip; i+=skip){
        ctxgraph.beginPath();
        ctxgraph.moveTo(x2px(xl[i], xmin, dx), y2px(yl[i], ymin, dy));
        ctxgraph.lineTo(x2px(xl[i+skip], xmin, dx), y2px(yl[i+skip], ymin, dy));
        ctxgraph.stroke();
    }
}

function change_graph(){
    graph_type = document.getElementById('changegraph').value;
    draw_graph();
}

function calculateFlipTable(temp){
    wolfp = 1 - Math.exp( -2./temp );
}

/*function calculateFlipTable(temp){
    gtable_de = [];
    gtable_doflip = [];
    gtable_flipprob = [];
    for (var i=0; i<5; i++){
        de = -2*(2*i - 4) - gfield;
        arg = -de / temp;
        gtable_de[i] = de;
        gtable_doflip[i] = 1*(de<=0);
        gtable_flipprob[i] = Math.exp(arg) * (temp > 0);
    }
    for (var i=0; i<5; i++){
        de = 2*(2*i - 4) + gfield;
        arg = -de / temp;
        gtable_de[i+5] = de;
        gtable_doflip[i+5] = 1*(de<=0);
        gtable_flipprob[i+5] = Math.exp(arg) * (temp > 0);
    }

    wolfp = 1 - Math.exp( -2./temp );
}

function update_metropolis(){
    var x = Math.floor(Math.random()*gN);
    var y = Math.floor(Math.random()*gN);
    var ind = x + y*gN;
    var neigh = neighborCount(x, y, gN, gboard);

    var ind2 = Math.round(neigh + 5*(gboard[ind]+1)/2);
    if (gtable_doflip[ind2] || Math.random() < gtable_flipprob[ind2]){
        gboard[ind] = -gboard[ind];
        put_pixel(x, y, gpx_size, gboard[x+y*gN]);

        genergy += 2.0*de/(gN*gN);
        gmag += 2.0*gboard[ind]/(gN*gN);
    }
    gt += 1.0/(gN*gN);
}*/

/*===============================================================================
    initialization and drawing 
================================================================================*/
function clear(){
    ctx.fillStyle = 'rgba(200,200,200,0.2)';
    ctx.clearRect(0, 0, c.width, c.height);
    ctx.fillRect(0,0,c.width,c.height);
}

function cleargraph(){
    ctxgraph.fillStyle = 'rgba(200,200,200,0.2)';
    ctxgraph.clearRect(0, 0, c2.width, c2.height);
    ctxgraph.fillRect(0,0,c2.width,c2.height);
}

var tick = function(T) {
    var skip = frameskip;
    if (update_func=='metropolis') {
      skip = skip*gN*gN;
    }
    if (dodraw == true) {
        for (var i=0; i<skip; i++){
            frame++;
            update();
        }
        sweeps = 1.0*frame / (gN*gN);
        draw_all();
        requestAnimationFrame(tick, c);
    }
};

function change_num(){
    gN = parseInt(document.getElementById('changenum').value);
    init_board(gN, null);
}



var init = function() {
    // create the canvas element
    empty = document.createElement('canvas');
    empty.width = empty.height = 1;
    c = document.getElementById('canvas');
    c.style.cursor = 'url('+empty.toDataURL()+')';
    ctx = c.getContext('2d');
    c2 = document.getElementById('canvas-graph');
    c2.style.cursor = 'url('+empty.toDataURL()+')';
    ctxgraph = c2.getContext('2d');
    gbuffer = ctx.getImageData(0, 0, canvasN, canvasN);
    gbufferdata = gbuffer.data;

    calculateFlipTable(gT);

    Number.prototype.mod = function(n) {
        return ((this%n)+n)%n;
    }

    document.getElementById('label_temp_input').addEventListener("keydown", function(e) {
        if (e.keyCode == 13){ 
            e.preventDefault();
            step = document.getElementById('temp').step;
            min = document.getElementById('temp').min;
            tval = parseFloat(document.getElementById('label_temp_input').value);

            if (tval <= Math.pow(10, min))
                logval = min - 2*step;
            else
                logval = log10(tval);

            document.getElementById('temp').value = logval;
            update_temp();
            undotextbox('label_temp_input');
        }
    }, false);

    document.getElementById('label_field_input').addEventListener("keydown", function(e) {
        if (e.keyCode == 13){ 
            e.preventDefault();
            document.getElementById('field').value = document.getElementById('label_field_input').value;
            update_field();
            undotextbox('label_field_input');
        }
    }, false);

    document.getElementById('label_frames_input').addEventListener("keydown", function(e) {
        if (e.keyCode == 13){ 
            e.preventDefault();
            tval = parseFloat(document.getElementById('label_frames_input').value);
            if (update_func=='metropolis')
                document.getElementById('frames').value = log10(tval);
            else
                document.getElementById('frames').value = tval;
            update_frames();
            undotextbox('label_frames_input');
        }
    }, false);

    clear();
    cleargraph();
    init_board(gN, null);
    update_display();

    document.body.addEventListener('keyup', function(ev) {
        if (ev.keyCode == 32){ ev.preventDefault(); update_pause(); } //space is pause
    }, false);

    document.body.addEventListener('keydown', function(ev) {
    }, false);

    registerAnimationRequest();
    requestAnimationFrame(tick, c);
};
window.onload = init;


// Provides requestAnimationFrame in a cross browser way.
// http://paulirish.com/2011/requestanimationframe-for-smart-animating/
function registerAnimationRequest() {
if ( !window.requestAnimationFrame ) {
    window.requestAnimationFrame = ( function() {
      return window.webkitRequestAnimationFrame ||
      window.mozRequestAnimationFrame || // comment out if FF4 is slow (it caps framerate at ~30fps)
      window.oRequestAnimationFrame ||
      window.msRequestAnimationFrame ||
      function( /* function FrameRequestCallback */ callback, /* DOMElement Element */ element ) {
              window.setTimeout( callback, 1 ); /*1000 / 60 );*/
      };
    } )();
}
}


