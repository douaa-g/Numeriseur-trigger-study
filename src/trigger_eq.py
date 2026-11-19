"""Trigger-equation compiler: text equation -> 16-bit truth table for four discriminated inputs A, B, C, D.

The design idea is the one used by the LPSC UCTM / UCTM2 modules (arXiv:1603.06816): the user types a logic
equation, software evaluates it for the 16 input combinations and loads the resulting truth table into a small
FPGA memory at run start, so the firmware is never resynthesised.

Grammar (case-insensitive), lowest to highest precedence:
    or_expr   := xor_expr  (( OR | NOR  | '|' ) xor_expr)*
    xor_expr  := and_expr  (( XOR | XNOR | '^' ) and_expr)*
    and_expr  := not_expr  (( AND | NAND | '&' ) not_expr)*
    not_expr  := ( NOT | '!' | '~' ) not_expr | atom
    atom      := A | B | C | D | '(' or_expr ')' | MULT '(' k ',' name {',' name} ')'      # at least k of the listed inputs
Truth-table bit index = A + 2*B + 4*C + 8*D (A is the least significant input).
Entry 0 (no input fired) is never used: the hardware gate only opens on a hit."""
import re

NAMES = "ABCD"
_TOK = re.compile(r"\s*(?:(?P<num>\d+)|(?P<id>[A-Za-z_]+)|(?P<sym>[()&|^!~,]))")

def _tokens(s):
    pos, out = 0, []
    s = s.strip()
    while pos < len(s):
        m = _TOK.match(s, pos)
        if not m:
            raise ValueError(f"unexpected character at position {pos}: {s[pos:pos+10]!r}")
        pos = m.end()
        if m.group("num"):
            out.append(("num", int(m.group("num"))))
        elif m.group("id"):
            out.append(("id", m.group("id").upper()))
        else:
            out.append(("sym", m.group("sym")))
    return out

class _P:
    def __init__(self, toks):
        self.t, self.i = toks, 0
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)
    def eat(self, kind=None, val=None):
        k, v = self.peek()
        if (kind and k != kind) or (val is not None and v != val):
            raise ValueError(f"expected {val or kind}, found {v}")
        self.i += 1
        return v
    def _bin(self, sub, ops):
        left = sub()
        while True:
            k, v = self.peek()
            tag = v if (k in ("id", "sym") and v in ops) else None
            if tag is None:
                return left
            self.i += 1
            right = sub()
            left = (ops[tag], left, right)
    def or_expr(self):
        return self._bin(self.xor_expr, {"OR": "or", "|": "or", "NOR": "nor"})
    def xor_expr(self):
        return self._bin(self.and_expr, {"XOR": "xor", "^": "xor", "XNOR": "xnor"})
    def and_expr(self):
        return self._bin(self.not_expr, {"AND": "and", "&": "and", "NAND": "nand"})
    def not_expr(self):
        k, v = self.peek()
        if v in ("NOT", "!", "~"):
            self.i += 1
            return ("not", self.not_expr())
        return self.atom()
    def atom(self):
        k, v = self.peek()
        if k == "id" and v in NAMES:
            self.i += 1
            return ("in", NAMES.index(v))
        if k == "id" and v == "MULT":
            self.i += 1
            self.eat("sym", "(")
            kmin = self.eat("num")
            idx = []
            while self.peek() == ("sym", ","):
                self.i += 1
                n = self.eat("id")
                if n not in NAMES:
                    raise ValueError(f"unknown input {n}")
                idx.append(NAMES.index(n))
            self.eat("sym", ")")
            if not idx or not 1 <= kmin <= len(idx):
                raise ValueError("MULT(k, inputs...) needs 1 <= k <= number of inputs")
            return ("mult", kmin, tuple(idx))
        if (k, v) == ("sym", "("):
            self.i += 1
            e = self.or_expr()
            self.eat("sym", ")")
            return e
        if k is None:
            raise ValueError("incomplete expression")
        raise ValueError(f"unexpected token {v!r}")

def _eval(node, bits):
    op = node[0]
    if op == "in":   return bits[node[1]]
    if op == "not":  return 1 - _eval(node[1], bits)
    if op == "mult": return int(sum(bits[i] for i in node[2]) >= node[1])
    a, b = _eval(node[1], bits), _eval(node[2], bits)
    return {"and": a & b, "or": a | b, "xor": a ^ b, "nand": 1 - (a & b),
            "nor": 1 - (a | b), "xnor": 1 - (a ^ b)}[op]

def compile_equation(eq):
    """Returns the 16-bit truth table (int). Raises ValueError on syntax errors."""
    p = _P(_tokens(eq))
    tree = p.or_expr()
    if p.peek()[0] is not None:
        raise ValueError(f"trailing input after expression: {p.peek()[1]!r}")
    tt = 0
    for idx in range(16):
        bits = [(idx >> i) & 1 for i in range(4)]
        tt |= _eval(tree, bits) << idx
    tt &= 0xFFFE            # entry 0 is meaningless (the gate opens on a hit), keep it clear
    return tt

def describe(tt):
    return [("".join(n for i, n in enumerate(NAMES) if (idx >> i) & 1)) for idx in range(16) if (tt >> idx) & 1]

if __name__ == "__main__":
    def table(f):
        return sum(1 << i for i in range(16) if f([(i >> k) & 1 for k in range(4)]))
    assert compile_equation("A & B") == 0x8888
    assert compile_equation("A AND B AND C") == 0x8080
    assert compile_equation("(A & B) | (C & D)") == 0xF888
    assert compile_equation("A XOR B") == 0x6666
    assert compile_equation("MULT(2, A, B, C, D)") == table(lambda b: sum(b) >= 2)
    assert compile_equation("A & B & !D") == table(lambda b: b[0] and b[1] and not b[3])
    assert compile_equation("a and b and not d") == compile_equation("A & B & !D")
    assert compile_equation("A OR B") == compile_equation("A | B")
    for bad in ["A & ", "E & A", "MULT(5, A, B)", "A B", "A &&& B"]:
        try:
            compile_equation(bad)
        except ValueError:
            continue
        raise SystemExit(f"should have been rejected: {bad}")
    assert compile_equation("B & !A") == table(lambda b: b[1] and not b[0])
    for eq in ["A & B", "A & B & C", "(A & B) & !D", "MULT(2, A, B, C, D)", "A XOR B"]:
        print(f"{eq:24s} -> 0x{compile_equation(eq):04X}   fires on: {describe(compile_equation(eq))}")
    print("trigger_eq self-test passed")
