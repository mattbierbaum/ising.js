// all of the global variables for dynamics
var gpx_black = null;
var gpx_white = null;
var gpx_size = 0;
var gboard = null;
var gN = 256;
var gt = 0;
var gT = 2.26918531421;
var gT = 2.28;
var gfield = 0;
var canvasN = 512;
var gbuffer;
var gbufferdata;

var gtable_doflip;
var gtable_flipprob;
var wolfp = 1-Math.exp(-2./gT);

// display variables
var c, c2;
var ctx;
var ctxgraph;
var empty;
var frame = 0;
var keys = [0,0,0,0];
var frameskip = 0.5;
var dodraw = true;

function rgb(r,g,b) {
    return 'rgb('+r+','+g+','+b+')';
}


function toFixed(value, precision) {
    var precision = precision || 0,
    neg = value < 0,
    power = Math.pow(10, precision),
    value = Math.round(value * power),
    integral = String((neg ? Math.ceil : Math.floor)(value / power)),
    fraction = String((neg ? -value : value) % power),
    padding = new Array(Math.max(precision - fraction.length, 0) + 1).join('0');
    var sneg = neg ? "-" : "";
    return sneg + (precision ? integral + '.' +  padding + fraction : integral);
}

function init_board(N, board){
    console.log("Beginning with "+N);
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

/*function energy(x, y, N, b){
    return 2*b[x+y*N]*(b[x + ((y+1).mod(N))*N] + 
        b[x + ((y-1).mod(N))*N] + 
        b[(x+1).mod(N) + y*N] + 
        b[(x-1).mod(N) + y*N] + gfield);
}


function update_metropolis(){
    var x = Math.floor(Math.random()*gN);
    var y = Math.floor(Math.random()*gN);
    var ind = x + y*gN;
    var de = energy(x, y, gN, gboard);
    if (de < 0 || Math.random() < Math.exp(-de / gT)){
        if (gboard[ind] == 1) 
            gboard[ind] = -1;
        else 
            gboard[ind] = 1;
        put_pixel(x, y, gpx_size, gboard[x+y*gN]);
    }
}*/

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
    var ds = 2 * state * Object.keys(cluster).length * gfield;

    if ( (ds < 0) || (Math.random() < Math.exp(-ds)) ) {
        // flip the cluster
        for (var ind in cluster) {
            ind = Number(ind);
            x = ind % gN;
            y = Math.floor( ind / gN );
            gboard[ind] = -state;
            put_pixel(x,y, gpx_size, -state);
        }
    }
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


function draw_all(){
    gbuffer.data = gbufferdata;
    ctx.putImageData(gbuffer, 0, 0);
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
    gT = parseFloat(document.getElementById('temp').value);
    document.getElementById('label_temp').innerHTML = toFixed(gT,6);
    calculateFlipTable(gT);
}
function update_field(){
    gfield = parseFloat(document.getElementById('field').value);
    document.getElementById('label_field').innerHTML = toFixed(gfield,6);
    calculateFlipTable(gT);
}
function update_frames(){
    frameval = parseFloat(document.getElementById('frames').value);
    if (update_func=='metropolis') {
        frameskip = Math.pow(10, frameval);
    } else {
      frameskip = frameval;
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
        frameskip = Math.pow(10.,-0.5);
        frame_label.innerHTML = toFixed(0.5,2);
        frame_slider.step = 0.01;
        frame_slider.max=0.5;
        frame_slider.min=-2;
        frame_slider.value = -0.5;
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

/*===============================================================================
 * graphing
 *=============================================================================*/
/*function calculateTics(xmax, xmin, ymax, ymin){
    var dx = xmax - xmin;
    var dy = ymax - ymin;

    var oom_x = Math.abs(dx)<1?Math.round(Math.log10(dx)):Math.floor(Math.log10(dx));
    var oom_y = Math.abs(dy)<1?Math.round(Math.log10(dy)):Math.floor(Math.log10(dy));
    var pow10_x = Math.pow(10, oom_x);
    var pow10_y = Math.pow(10, oom_y);

    int idx = (int)(dx / pow10_x);
    int idy = (int)(dy / pow10_y);

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
}*/
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
    var neigh = neighborCount(x, y, gN, gboard);

    var ind2 = Math.round(neigh + 5*(gboard[ind]+1)/2);
    if (gtable_doflip[ind2] || Math.random() < gtable_flipprob[ind2]){
        if (gboard[ind] == 1) 
            gboard[ind] = -1;
        else 
            gboard[ind] = 1;
        put_pixel(x, y, gpx_size, gboard[x+y*gN]);
    }
}

/*===============================================================================
    initialization and drawing 
================================================================================*/
function clear(){
    ctx.fillStyle = 'rgba(200,200,200,0.2)';
    ctx.clearRect(0, 0, c.width, c.height);
    ctx.fillRect(0,0,c.width,c.height);
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
        draw_all();
        requestAnimationFrame(tick, c);
    }
};

function change_num(){
    gN = parseInt(document.getElementById('changenum').value);
    init_board(gN, null);
}

function calculateFlipTable(temp){
    gtable_doflip = [];
    gtable_flipprob = [];
    for (var i=0; i<5; i++){
        de = -2*(2*i - 4) - gfield;
        arg = -de / temp;
        gtable_doflip[i] = 1*(de<=0);
        gtable_flipprob[i] = Math.exp(arg) * (temp > 0);
    }
    for (var i=0; i<5; i++){
        de = 2*(2*i - 4) + gfield;
        arg = -de / temp;
        gtable_doflip[i+5] = 1*(de<=0);
        gtable_flipprob[i+5] = Math.exp(arg) * (temp > 0);
    }

    wolfp = 1 - Math.exp( -2./temp );
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
            document.getElementById('temp').value = document.getElementById('label_temp_input').value;
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
            document.getElementById('label_frames').innerHTML = document.getElementById('label_frames_input').value;
            update_frames();
            undotextbox('label_frames_input');
        }
    }, false);

    clear();
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


