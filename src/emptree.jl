"""
Add a node of type MP


# Arguments
- emp: the EMP object
- jamsd_emp: an opaque pointer to EMP object for jamsd
- node: the node to process
- ctx: an opaque pointer to the context object for jamsd
- jamsd_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMP.Model, jamsd_emp::Ptr{JAMSDWriter.empinfo}, node::MathPrgm, ctx::Ptr{JAMSDWriter.context}, jamsd_mp_parent::Ptr{JAMSDWriter.mathprgm}=Ptr{JAMSDWriter.mathprgm}(C_NULL))
    jamsd_mp = JAMSDWriter.jamsd_declare_mathprgm(node, ctx, jamsd_emp)
    node.solverobj = jamsd_mp

#    JAMSDWriter.jamsd_add_mp(jamsd_emp, jamsd_node, jamsd_mp)

    if jamsd_mp_parent != C_NULL
        JAMSDWriter.emp_add_mp_mp(jamsd_mp_parent, jamsd_mp)
    else
        JAMSDWriter.emp_set_root(jamsd_emp, jamsd_mp)
    end

    for mp_mp in node.mps
        empnode(emp, jamsd_emp, mp_mp, ctx, jamsd_mp)
    end

    for mp_equil in node.equils
        empnode(emp, jamsd_emp, mp_equil, ctx, jamsd_mp)
    end

end

"""
Add a node of type equilibrium


# Arguments
- emp: the EMP object
- jamsd_emp: an opaque pointer to EMP object for jamsd
- node: the node to process
- ctx: an opaque pointer to the context object for jamsd
- jamsd_mp_parent: an opaque pointer to the parent of the node
"""
function empnode(emp::EMP.Model, jamsd_emp::Ptr{JAMSDWriter.empinfo}, node::Vector{MathPrgm}, ctx::Ptr{JAMSDWriter.context}, jamsd_mp_parent::Ptr{JAMSDWriter.mathprgm}=Ptr{JAMSDWriter.mathprgm}(C_NULL))
    jamsd_mpe = JAMSDWriter.emp_create_equil(length(node))

    jamsd_mps = []
    for mp in node
        jamsd_mp = JAMSDWriter.jamsd_declare_mathprgm(mp, ctx, jamsd_emp)
        mp.solverobj = jamsd_mp
        JAMSDWriter.emp_equil_add(jamsd_mpe, jamsd_mp)
        push!(jamsd_mps, jamsd_mp)
    end

    if jamsd_mp_parent != C_NULL
        JAMSDWriter.emp_add_mp_equil(jamsd_emp, jamsd_node, jamsd_mpe)
    else
        JAMSDWriter.emp_set_root(jamsd_emp, jamsd_mpe)
    end

    for (mp, jamsd_mp) in zip(node, jamsd_mps)
        for mp_mp in mp.mps
            empnode(emp, jamsd_emp, mp_mp, ctx, jamsd_mp)
        end

        for mp_equil in mp.equils
            empnode(emp, jamsd_emp, mp_equil, ctx, jamsd_mp)
        end
    end
end

"""
    emptree(emp, jamsd_emp, ctx)

Pass the EMP graph structure to the EMP solver
"""
function emptree(emp::EMP.Model, jamsd_emp::Ptr{JAMSDWriter.empinfo}, ctx::Ptr{JAMSDWriter.context})
    if emp.root.hasvalue
        JAMSDWriter.emp_mp_ensure(jamsd_emp, length(emp.mps))

        empnode(emp, jamsd_emp, emp.root.value, ctx)
    end
end
