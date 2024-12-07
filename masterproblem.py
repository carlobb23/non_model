import gurobipy as gu

class MasterProblem:
    def __init__(self, df, max_iteration, current_iteration, T_Max):
        self.iteration = current_iteration
        self.P = df['P'].dropna().astype(int).unique().tolist()
        self.D = df['D'].dropna().astype(int).unique().tolist()
        self.T = df['T'].dropna().astype(int).unique().tolist()
        self.A = [i for i in range(1, max_iteration + 2)]
        self.Model = gu.Model("MasterProblem")
        self.cons_p_max = {}
        self.cons_lmbda = {}
        self.T_max = T_Max

    def buildModel(self):
        self.genVars()
        self.genCons()
        self.genObj()
        self.Model.update()

    def genVars(self):
        self.lmbda = self.Model.addVars(self.P, self.A, vtype=gu.GRB.BINARY,  name='lmbda')

    def genCons(self):
        for p in self.P:
            self.cons_lmbda[p] = self.Model.addConstr(1 == gu.quicksum(self.lmbda[p, a] for a in self.A),name=f"lambda({p})")
        for t in self.T:
            for d in self.D:
                self.cons_p_max[t, d] = self.Model.addConstr(
                    gu.quicksum(self.lmbda[p, a] for p in self.P for a in self.A) <= self.T_max[t, d], name=f"p_max({t},{d})")

    def genObj(self):
        self.Model.setObjective(gu.quicksum(self.lmbda[p, a] for p in self.P for a in self.A), sense=gu.GRB.MINIMIZE)
    def getDuals(self):
        return {(t, d): self.cons_p_max[t, d].Pi for t in self.T for d in self.D}, {p: self.cons_lmbda[p].Pi for p in self.P}

    def initCoeffs(self):
        for p in self.P:
            for a in self.A[1:]:
                self.lmbda[p, a].Obj = 1000
        for t in self.T:
            for d in self.D:
                for p in self.P:
                    for a in self.A[1:]:
                        self.Model.chgCoeff(self.cons_p_max[t, d], self.lmbda[p, a], 1000)
        self.Model.update()

    def startSol(self, schedules_x, schedules_los):
        for p in self.P:
            for t in self.T:
                for d in self.D:
                    if (p, t, d) in schedules_x:
                        value = schedules_x[p, t, d]
                    else:
                        value = 0
                    self.lmbda[p, 1].Obj = schedules_los.get((p), 0)

                    if (p, t, d) in schedules_x:
                        self.Model.chgCoeff(self.cons_p_max[t, d], self.lmbda[p, 1], value)
        self.Model.update()

    def addCol(self, p, it, schedules_x, schedules_o):
        iter = it + 1
        for t in self.T:
            for d in self.D:
                if (p, t, d, iter) in schedules_x:
                    value = schedules_x[p, t, d, iter]
                else:
                    value = 0
                    print(f'No!')

                self.lmbda[p, iter].Obj = schedules_o.get((p, iter), 0)

                if (p, t, d, iter) in schedules_x:
                    self.Model.chgCoeff(self.cons_p_max[t, d], self.lmbda[p, iter], value)
        self.Model.update()

    def finSol(self):
        self.Model.Params.OutputFlag = 0
        self.Model.setAttr("vType", self.lmbda, gu.GRB.INTEGER)
        self.Model.update()
        self.Model.optimize()

    def solRelModel(self):
        self.Model.Params.OutputFlag = 0
        for v in self.Model.getVars():
            v.setAttr('vtype', 'C')
            v.setAttr('lb', 0.0)
        self.Model.update()
        self.Model.optimize()