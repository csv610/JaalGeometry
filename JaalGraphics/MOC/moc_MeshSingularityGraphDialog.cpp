/****************************************************************************
** Meta object code from reading C++ file 'MeshSingularityGraphDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshSingularityGraphDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshSingularityGraphDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshSingularityGraphDialog_t {
    QByteArrayData data[10];
    char stringdata0[145];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshSingularityGraphDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshSingularityGraphDialog_t qt_meta_stringdata_JMeshSingularityGraphDialog = {
    {
QT_MOC_LITERAL(0, 0, 27), // "JMeshSingularityGraphDialog"
QT_MOC_LITERAL(1, 28, 24), // "openNodeAttributesDialog"
QT_MOC_LITERAL(2, 53, 0), // ""
QT_MOC_LITERAL(3, 54, 24), // "openEdgeAttributesDialog"
QT_MOC_LITERAL(4, 79, 17), // "mouseReleaseEvent"
QT_MOC_LITERAL(5, 97, 12), // "QMouseEvent*"
QT_MOC_LITERAL(6, 110, 1), // "e"
QT_MOC_LITERAL(7, 112, 8), // "getGraph"
QT_MOC_LITERAL(8, 121, 11), // "displayMesh"
QT_MOC_LITERAL(9, 133, 11) // "closeDialog"

    },
    "JMeshSingularityGraphDialog\0"
    "openNodeAttributesDialog\0\0"
    "openEdgeAttributesDialog\0mouseReleaseEvent\0"
    "QMouseEvent*\0e\0getGraph\0displayMesh\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshSingularityGraphDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    1,   46,    2, 0x08 /* Private */,
       7,    0,   49,    2, 0x08 /* Private */,
       8,    0,   50,    2, 0x08 /* Private */,
       9,    0,   51,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    6,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshSingularityGraphDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshSingularityGraphDialog *_t = static_cast<JMeshSingularityGraphDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openNodeAttributesDialog(); break;
        case 1: _t->openEdgeAttributesDialog(); break;
        case 2: _t->mouseReleaseEvent((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 3: _t->getGraph(); break;
        case 4: _t->displayMesh(); break;
        case 5: _t->closeDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JMeshSingularityGraphDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshSingularityGraphDialog.data,
      qt_meta_data_JMeshSingularityGraphDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshSingularityGraphDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshSingularityGraphDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshSingularityGraphDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshSingularityGraphDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshSingularityGraphDialog"))
        return static_cast< Ui::MeshSingularityGraphDialog*>(const_cast< JMeshSingularityGraphDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshSingularityGraphDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
