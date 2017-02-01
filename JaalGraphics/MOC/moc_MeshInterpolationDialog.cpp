/****************************************************************************
** Meta object code from reading C++ file 'MeshInterpolationDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshInterpolationDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshInterpolationDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshInterpolationDialog_t {
    QByteArrayData data[11];
    char stringdata0[152];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshInterpolationDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshInterpolationDialog_t qt_meta_stringdata_JMeshInterpolationDialog = {
    {
        QT_MOC_LITERAL(0, 0, 24), // "JMeshInterpolationDialog"
        QT_MOC_LITERAL(1, 25, 12), // "checkDisplay"
        QT_MOC_LITERAL(2, 38, 0), // ""
        QT_MOC_LITERAL(3, 39, 12), // "getImageMesh"
        QT_MOC_LITERAL(4, 52, 13), // "openLIMDialog"
        QT_MOC_LITERAL(5, 66, 25), // "openMeshGeomQualityDialog"
        QT_MOC_LITERAL(6, 92, 11), // "closeDialog"
        QT_MOC_LITERAL(7, 104, 10), // "loadSource"
        QT_MOC_LITERAL(8, 115, 10), // "loadTarget"
        QT_MOC_LITERAL(9, 126, 12), // "deformSolver"
        QT_MOC_LITERAL(10, 139, 12) // "getHausdorff"

    },
    "JMeshInterpolationDialog\0checkDisplay\0"
    "\0getImageMesh\0openLIMDialog\0"
    "openMeshGeomQualityDialog\0closeDialog\0"
    "loadSource\0loadTarget\0deformSolver\0"
    "getHausdorff"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshInterpolationDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    9,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   59,    2, 0x08 /* Private */,
    3,    0,   60,    2, 0x08 /* Private */,
    4,    0,   61,    2, 0x08 /* Private */,
    5,    0,   62,    2, 0x08 /* Private */,
    6,    0,   63,    2, 0x08 /* Private */,
    7,    0,   64,    2, 0x08 /* Private */,
    8,    0,   65,    2, 0x08 /* Private */,
    9,    0,   66,    2, 0x08 /* Private */,
    10,    0,   67,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMeshInterpolationDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshInterpolationDialog *_t = static_cast<JMeshInterpolationDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->checkDisplay();
            break;
        case 1:
            _t->getImageMesh();
            break;
        case 2:
            _t->openLIMDialog();
            break;
        case 3:
            _t->openMeshGeomQualityDialog();
            break;
        case 4:
            _t->closeDialog();
            break;
        case 5:
            _t->loadSource();
            break;
        case 6:
            _t->loadTarget();
            break;
        case 7:
            _t->deformSolver();
            break;
        case 8:
            _t->getHausdorff();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshInterpolationDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshInterpolationDialog.data,
        qt_meta_data_JMeshInterpolationDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshInterpolationDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshInterpolationDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshInterpolationDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshInterpolationDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshInterpolationDialog"))
        return static_cast< Ui::MeshInterpolationDialog*>(const_cast< JMeshInterpolationDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshInterpolationDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
