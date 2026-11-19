// Digital constant-fraction discriminator (CFD), fixed point.
//   y[n] = (K*x[n] >>> 3) - x[n-D]        (fraction f = K/8, delay D samples)
// Arms when x >= THR, then reports the first positive-to-negative zero crossing of y
// (sample index n, and y[n-1], y[n]).  The sub-sample position is  (n-1) + y[n-1]/(y[n-1]-y[n])
// and is computed downstream (software / a small divider); the RTL stays divider-free.
// After a hit the unit is locked out until x has stayed below THR for HOLD consecutive samples
// (prevents re-triggering on noise in the pulse tail).
module cfd #(
    parameter integer W   = 12,   // ADC width (signed)
    parameter integer D   = 6,    // delay in samples
    parameter integer K   = 3,    // fraction K/8
    parameter integer THR = 90,   // arming threshold in LSB
    parameter integer HOLD = 16   // re-arm only after x < THR for HOLD consecutive samples
) (
    input  wire                   clk,
    input  wire                   rst,
    input  wire signed [W-1:0]    x,
    output reg                    valid,
    output reg  [31:0]            n_out,
    output reg  signed [W+3:0]    yp_out,   // y[n-1]
    output reg  signed [W+3:0]    yc_out    // y[n]
);
    reg signed [W-1:0] dl [0:D-1];          // delay line: dl[D-1] = x[n-D]
    reg [31:0]         n;
    reg                armed, lock;
    reg [7:0]          low_cnt;
    reg signed [W+3:0] y_prev;
    integer            i;

    wire signed [W+3:0] kx   = K * x;
    wire signed [W+3:0] kx_s = kx >>> 3;
    wire signed [W+3:0] xd   = dl[D-1];
    wire signed [W+3:0] y    = kx_s - xd;
    wire                fire = armed && (y_prev > 0) && (y <= 0);

    always @(posedge clk) begin
        if (rst) begin
            n <= 0; armed <= 0; lock <= 0; low_cnt <= 0; y_prev <= 0;
            valid <= 0; n_out <= 0; yp_out <= 0; yc_out <= 0;
            for (i = 0; i < D; i = i + 1) dl[i] <= 0;
        end else begin
            n     <= n + 1;
            valid <= fire;
            if (fire) begin
                n_out  <= n;
                yp_out <= y_prev;
                yc_out <= y;
            end
            if (lock) begin
                if (x < THR) begin
                    if (low_cnt == HOLD - 1) begin lock <= 0; low_cnt <= 0; end
                    else low_cnt <= low_cnt + 1;
                end else low_cnt <= 0;
            end else if (armed) begin
                if (fire) begin armed <= 0; lock <= 1; end
            end else if (x >= THR) begin
                armed <= 1;
            end
            y_prev <= y;
            dl[0]  <= x;
            for (i = 1; i < D; i = i + 1) dl[i] <= dl[i-1];
        end
    end
endmodule
