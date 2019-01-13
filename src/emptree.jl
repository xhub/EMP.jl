"""
Add a node of type MP


# Arguments
- emp: the EMP object
- reshop_emp: an opaque pointer to EMP object for ReSHOP
- node: the node to process
- ctx: an opaque pointer to the context object for ReSHOP
- reshop_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMP.Model, reshop_emp::Ptr{ReSHOP.empinfo}, node::MathPrgm, ctx::Ptr{ReSHOP.context}, reshop_mp_parent::Ptr{ReSHOP.mathprgm}=Ptr{ReSHOP.mathprgm}(C_NULL))
    reshop_mp = ReSHOP.reshop_declare_mathprgm(node, ctx, reshop_emp)
    node.solverobj = reshop_mp

#    ReSHOP.reshop_add_mp(reshop_emp, reshop_node, reshop_mp)

    if reshop_mp_parent != C_NULL
        ReSHOP.emp_add_mp_mp(reshop_mp_parent, reshop_mp)
    else
        ReSHOP.emp_set_root(reshop_emp, reshop_mp)
    end

    for mp_mp in node.mps
        empnode(emp, reshop_emp, mp_mp, ctx, reshop_mp)
    end

    for mp_equil in node.equils
        empnode(emp, reshop_emp, mp_equil, ctx, reshop_mp)
    end

end

"""
Add a node of type equilibrium


# Arguments
- emp: the EMP object
- reshop_emp: an opaque pointer to EMP object for ReSHOP
- node: the node to process
- ctx: an opaque pointer to the context object for ReSHOP
- reshop_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMP.Model, reshop_emp::Ptr{ReSHOP.empinfo}, node::Vector{MathPrgm}, ctx::Ptr{ReSHOP.context}, reshop_mp_parent::Ptr{ReSHOP.mathprgm}=Ptr{ReSHOP.mathprgm}(C_NULL))
    reshop_mpe = ReSHOP.emp_create_equil(length(node))

    reshop_mps = []
    for mp in node
        reshop_mp = ReSHOP.reshop_declare_mathprgm(mp, ctx, reshop_emp)
        mp.solverobj = reshop_mp
        ReSHOP.emp_equil_add(reshop_mpe, reshop_mp)
        push!(reshop_mps, reshop_mp)
    end

    if reshop_mp_parent != C_NULL
        ReSHOP.emp_add_mp_equil(reshop_emp, reshop_node, reshop_mpe)
    else
        ReSHOP.emp_set_root(reshop_emp, reshop_mpe)
    end

    for (mp, reshop_mp) in zip(node, reshop_mps)
        for mp_mp in mp.mps
            empnode(emp, reshop_emp, mp_mp, ctx, reshop_mp)
        end

        for mp_equil in mp.equils
            empnode(emp, reshop_emp, mp_equil, ctx, reshop_mp)
        end
    end
end

"""
    emptree(emp, reshop_emp, ctx)

Pass the EMP graph structure to the EMP solver
"""
function emptree(emp::EMP.Model, reshop_emp::Ptr{ReSHOP.empinfo}, ctx::Ptr{ReSHOP.context})
    if emp.root.hasvalue
        ReSHOP.emp_mp_ensure(reshop_emp, length(emp.mps))

        empnode(emp, reshop_emp, emp.root.value, ctx)
    end
end
