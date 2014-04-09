// all of the global variables for dynamics
var gpx_black = null;
var gpx_white = null;
var gpx_size = 0;
var gboard = null;
var gN = 400;
var gt = 0;
var gT = 2.26918531421;
var gfield = 0;
var canvasN = 400;
var gbuffer;
var gbufferdata;

// display variables
var c;
var ctx;
var empty;
var frame = 0;
var keys = [0,0,0,0];
var frameskip = 5000;
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

function energy(x, y, N, b){
    return 2*b[x+y*N]*(b[x + ((y+1).mod(N))*N] + 
        b[x + ((y-1).mod(N))*N] + 
        b[(x+1).mod(N) + y*N] + 
        b[(x-1).mod(N) + y*N] + gfield);
}

function update(){
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
}

function draw_all(){
    gbuffer.data = gbufferdata;
    ctx.putImageData(gbuffer, 0, 0);
}


/*======================================================================
  the javascript interface stuff
=========================================================================*/
function update_temp(){
    gT = parseFloat(document.getElementById('temp').value);
    document.getElementById('label_temp').innerHTML = toFixed(gT,2);
}
function update_field(){
    gfield = parseFloat(document.getElementById('field').value);
    document.getElementById('label_field').innerHTML = toFixed(gfield,2);
}
function update_frames(){
    frameskip = document.getElementById('frames').value;
    document.getElementById('label_frames').innerHTML = toFixed(frameskip,2);
}

function update_display(){
    document.getElementById('label_temp').innerHTML = toFixed(gT,2);
    document.getElementById('label_field').innerHTML = toFixed(gfield,2);
    document.getElementById('label_frames').innerHTML = toFixed(frameskip,2);
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
    initialization and drawing 
================================================================================*/
function clear(){
    ctx.fillStyle = 'rgba(200,200,200,0.2)';
    ctx.clearRect(0, 0, c.width, c.height);
    ctx.fillRect(0,0,c.width,c.height);
}

var tick = function(T) {
    if (dodraw == true) {
        for (var i=0; i<frameskip; i++){
            frame++;
            update();
        }
        draw_all();
        requestAnimationFrame(tick, c);
    }
};

function update_allcontrols(){
    document.getElementById('num').value = gN;
    document.getElementById('num').innerHTML = toFixed(0,2);
}

function change_num(){
    gN = parseInt(document.getElementById('num').value);
    init_board(gN, null);
}

var init = function() {
    // create the canvas element
    empty = document.createElement('canvas');
    empty.width = empty.height = 1;
    c = document.getElementById('canvas');
    c.style.cursor = 'url('+empty.toDataURL()+')';
    ctx = c.getContext('2d');
    gbuffer = ctx.getImageData(0, 0, canvasN, canvasN);
    gbufferdata = gbuffer.data;

    Number.prototype.mod = function(n) {
        return ((this%n)+n)%n;
    }

    clear();
    init_board(400, null);
    update_display();
    update_allcontrols();

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


